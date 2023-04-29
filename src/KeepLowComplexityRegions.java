import java.util.Arrays;
import java.io.*;


/**
 * Prints to output all maximal low-complexity intervals of the chromosomes.
 * Join two calls by an edge if they overlap. A connected component of this
 * graph is "complex" iff there are two starting positions or two ending 
 * positions in it that are farther than a threshold. The projection of a
 * complex component on its chromosome is marked as a complex interval.
 */
public class KeepLowComplexityRegions {
    /**
     * Total lengths
     */
    private static final int CHR21_LENGTH = 46709983;
    private static final int CHR22_LENGTH = 50818468;
    
    /**
     * SV types. INS are represented as intervals of constant length centered at
     * the insertion position.
     */
    private static final String[] SV_TYPES = new String[] {"del","inv","dup","ins","all"};
    private static final int N_SV_TYPES = 5;
    private static final int INS_SLACK = 100;  // Arbitrary
    
    /**
     * Calls
     */
    private static Call[] calls;
    private static int lastCall;
    
    /**
     * Connected components
     */
    private static int[] componentChr;
    private static int[] componentType, componentSize;
    private static int lastComponent;
    private static int[][] endpointMin, endpointMax;
    private static int[] lengthMin, lengthMax;
    
    
    /**
     * Remark: only INS and DEL are used to define complex intervals.
     *
     * @param args 2: consider as unreliable any DEL that is longer than this.
     */
    public static void main(String[] args) throws Exception {
        final String VCF_FILE = args[0];
        final boolean ONLY_PASS = Integer.parseInt(args[1])==1;
        final int MAX_DEL_LENGTH = Integer.parseInt(args[2]);
        final int DISTANCE_THRESHOLD = Integer.parseInt(args[3]);
        final int MIN_SIZE_COMPLEX = Integer.parseInt(args[4]);
        final String OUTPUT_DIR = args[5];
        
        System.err.print("Loading all calls... ");
        loadCalls(VCF_FILE,ONLY_PASS,MAX_DEL_LENGTH);
        System.err.println("DONE. "+(lastCall+1)+" calls loaded.");
        System.err.print("Building connected components... ");
        buildComponents();
        System.err.println("DONE. "+(lastComponent+1)+" components.");
        printIntervals(DISTANCE_THRESHOLD,MIN_SIZE_COMPLEX,OUTPUT_DIR+"/intervals_keep.txt");
    }
    
    
    /**
     * Remark: only INS and DEL are loaded.
     *
     * Remark: INS are represented here as intervals [pos..pos+length-1] for 
     * taking into account lengths later for deciding if a component is complex.
     *
     * Remark: the elements of $calls$ are not necessarily sorted when the
     * procedure completes.
     *
     * @param maxDelLength consider as unreliable any DEL longer than this.
     */
    private static final void loadCalls(String path, boolean onlyPass, int maxDelLength) throws IOException {
        final int CAPACITY = 1000;  // Arbitrary
        int i, n;
        int position, length, row, to, chr;
        String str;
        BufferedReader br;
        String[] tokens;
        
        calls = new Call[CAPACITY]; lastCall=-1;
        br = new BufferedReader(new FileReader(path));
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)==COMMENT) {
                str=br.readLine();
                continue;
            }
            tokens=str.split("\t");
            if (!tokens[0].equalsIgnoreCase("chr21") && !tokens[0].equalsIgnoreCase("chr22")) {
				str=br.readLine();
				continue;
            }
            position=Integer.parseInt(tokens[1]);
			if (onlyPass && !tokens[6].equalsIgnoreCase(PASS_STR)) {
				str=br.readLine();
				continue;
			}
            row=svType2Row(getField(tokens[7],SVTYPE_STR));
            if (row==-1) {
				str=br.readLine();
				continue;
            }
            length=Integer.parseInt(getField(tokens[7],SVLEN_STR));
            if (length<0) length=-length;
            if (row==0 && length>maxDelLength) {
                str=br.readLine();
                continue;
            }
            to=position+length-1;  // Yes, also for INS.
            if (tokens[0].equals("chr21")) chr=21;
            else if (tokens[0].equals("chr22")) chr=22;
            else chr=-1;
            lastCall++;
            if (lastCall==calls.length) {
                Call[] newArray = new Call[calls.length<<1];
                System.arraycopy(calls,0,newArray,0,calls.length);
                calls=newArray;
            }
            calls[lastCall] = new Call(row,chr,position,to);
            str=br.readLine();
        }
        br.close();
    }
    
    
    /**
     * Computes the connected component of every call in $calls$.
     *
     * Remark: two INS are connected iff the fixed-size intervals around their
     * starting positions overlap.
     *
     * Remark: assume that there is a simple cluster of long deletions, and 
     * assume that it covers several simple clusters of short deletions. All 
     * such long and short deletions would be put in the same component, which
     * would be marked as complex even though the actual clusters are simple.
     * Another similarly undesirable case occurs when a simple cluster covers a
     * complex cluster (e.g. because it covers a tandem).
     */
    private static final void buildComponents() {
        int i, j;
        int chr, type, first, last, component;
        
        System.err.print("Sorting calls...");
        Arrays.sort(calls,0,lastCall+1);
        System.err.println("DONE");
        lastComponent=-1;
        for (i=0; i<=lastCall; i++) {
            if (calls[i].component==-1) calls[i].component=++lastComponent;
            chr=calls[i].chr; type=calls[i].type;
            if (type==3) last=calls[i].first+INS_SLACK;
            else last=calls[i].last;
            for (j=i+1; j<=lastCall; j++) {
                if (calls[j].chr!=chr) break;
                if (calls[j].type!=type) continue;
                if (type==3) first=calls[j].first-INS_SLACK;
                else first=calls[j].first;
                if (first>last) break;
                calls[j].component=calls[i].component;
            }
        }
        componentChr = new int[lastComponent+1];
        componentType = new int[lastComponent+1];
        componentSize = new int[lastComponent+1];
        Arrays.fill(componentSize,0);
        for (i=0; i<=lastCall; i++) {
            componentChr[calls[i].component]=calls[i].chr;
            componentType[calls[i].component]=calls[i].type;
            componentSize[calls[i].component]++;
        }
        endpointMin = new int[lastComponent+1][2];
        for (i=0; i<=lastComponent; i++) { endpointMin[i][0]=Integer.MAX_VALUE; endpointMin[i][1]=Integer.MAX_VALUE; }
        endpointMax = new int[lastComponent+1][2];
        for (i=0; i<=lastCall; i++) {
            component=calls[i].component;
            if (calls[i].first<endpointMin[component][0]) endpointMin[component][0]=calls[i].first;
            if (calls[i].first>endpointMax[component][0]) endpointMax[component][0]=calls[i].first;
            if (calls[i].last<endpointMin[component][1]) endpointMin[component][1]=calls[i].last;
            if (calls[i].last>endpointMax[component][1]) endpointMax[component][1]=calls[i].last;
        }
    }
    
    
    /**
     * Writes to $path$ every maximal interval that is not covered by a complex
     * component. A component is complex iff it contains $>=minSizeComplex$
     * calls and it has a pair of start positions or end positions at distance
     * $>threshold$.
     *
     * Remark: for INS this corresponds to comparing start position and length;
     * the masking, however, is perfomed only around the start position.
     */
    private static final void printIntervals(int threshold, int minSizeComplex, String path) throws IOException {
        int i, j;
        int first, last;
        BufferedWriter bw;
        boolean[] mask, mask_chr21, mask_chr22;
        
        // Masking regions of the genome covered by some complex component
        mask_chr21 = new boolean[CHR21_LENGTH];
        mask_chr22 = new boolean[CHR22_LENGTH];
        for (i=0; i<=lastComponent; i++) {
            if ( componentSize[i]>=minSizeComplex &&
                 ( endpointMax[i][0]-endpointMin[i][0]>threshold || 
                   endpointMax[i][1]-endpointMin[i][1]>threshold
                 )
               ) {
                if (componentChr[i]==21) mask=mask_chr21;
                else if (componentChr[i]==22) mask=mask_chr22;
                else mask=null;
                if (componentType[i]==3) { 
                    first=endpointMin[i][0]-INS_SLACK;  
                    last=endpointMax[i][0]+INS_SLACK;  
                }
                else {
                    first=endpointMin[i][0];
                    last=endpointMax[i][1];
                }
                for (j=first; j<=last; j++) mask[j]=true;
            }
        }
        
        // Printing maximal unmasked intervals
        bw = new BufferedWriter(new FileWriter(path));
        first=-1; last=-1;
        for (i=0; i<CHR21_LENGTH; i++) {
            if (!mask_chr21[i]) {
                if (first==-1) first=i;
                last=i;
            }
            else {
                if (first!=-1) bw.write("chr21,"+first+","+last+"\n");
                first=-1; last=-1;
            }
        }
        if (first!=-1) bw.write("chr21,"+first+","+last+"\n");
        first=-1; last=-1;
        for (i=0; i<CHR22_LENGTH; i++) {
            if (!mask_chr22[i]) {
                if (first==-1) first=i;
                last=i;
            }
            else {
                if (first!=-1) bw.write("chr22,"+first+","+last+"\n");
                first=-1; last=-1;
            }
        }
        if (first!=-1) bw.write("chr22,"+first+","+last+"\n");
        bw.close();
    }
    
    
    private static class Call implements Comparable {
        public int type, chr, first, last, length;
        public int component;
        
        public Call() { }
        
        public Call(int t, int c, int f, int l) {
            this.type=t; this.chr=c; this.first=f; this.last=l;
            component=-1;
        }
        
        public int compareTo(Object other) {
            Call otherCall = (Call)other;
            if (chr<otherCall.chr) return -1;
            else if (chr>otherCall.chr) return 1;
            if (first<otherCall.first) return -1;
            else if (first>otherCall.first) return 1;
            if (last<otherCall.last) return -1;
            else if (last>otherCall.last) return 1;
            return 0;
        }
        
        public String toString() { 
            return chr+"\t"+first+"\t"+SV_TYPES[type]+"\t"+last+"\t"; 
        }
    }
    
    
    /**
     * This procedure returns just DEL and INS, since using the other types to
     * compute complex regions is not reliable with sniffles2.
     */
	private static final int svType2Row(String type) {
		if (type==null || type.length()==0) return -1;
		if ( type.equalsIgnoreCase(DEL_STR) || 
			 type.equalsIgnoreCase(DEL_ME_STR)
		   ) return 0;
		else if (type.equalsIgnoreCase(INV_STR)) return -1;
        else if ( type.equalsIgnoreCase(DUP_STR) ||
			      type.equalsIgnoreCase(DUP_TANDEM_STR) ||
				  type.equalsIgnoreCase(DUP_INT_STR)
			    ) return -1;
        else if ( type.equalsIgnoreCase(INS_STR) ||
                  type.equalsIgnoreCase(INS_ME_STR) ||
                  type.equalsIgnoreCase(INS_NOVEL_STR)
                ) return 3;
		else return -1;
	}
    
    
	/**
	 * @return NULL if $field$ does not occur in $str$.
	 */
	private static final String getField(String str, String field) {
		final int FIELD_LENGTH = field.length()+1;
		int p = str.indexOf(field+"=");
		if (p<0) return null;
		if (field.equalsIgnoreCase(END_STR)) {
			while (p>=2 && str.substring(p-2,p-2+CIEND_STR.length()).equalsIgnoreCase(CIEND_STR)) p=str.indexOf(field+"=",p+1);
			if (p<0) return null;
		}
		final int q = str.indexOf(SEPARATOR,p+FIELD_LENGTH);
		return str.substring(p+FIELD_LENGTH,q<0?str.length():q);
	}    
    
    
	/**
	 * Basic constants
	 */
	public static final char COMMENT = '#';
	public static final String SEPARATOR = ";";
	public static final String PASS_STR = "PASS";
	public static final String PRECISE_STR = "PRECISE";
	public static final String IMPRECISE_STR = "IMPRECISE";
	public static final String END_STR = "END";
	public static final String SVTYPE_STR = "SVTYPE";
	public static final String SVLEN_STR = "SVLEN";
	public static final String CHR2_STR = "CHR2";
	public static final String CT_STR = "CT";
	public static final String CT_325_STR = "'3to5'";
	public static final String CT_523_STR = "'5to3'";
	public static final String CT_525_STR = "'5to5'";
	public static final String CT_323_STR = "'3to3'";
    
	/**
	 * SV types: labels used by callers.
	 */
	public static final String DEL_STR = "DEL";
	public static final String DEL_ME_STR = "DEL:ME";
	public static final String DEL_INV_STR = "DEL/INV";
	public static final String INS_STR = "INS";
	public static final String INS_ME_STR = "INS:ME";
	public static final String INS_NOVEL_STR = "INS:NOVEL";
	public static final String DUP_STR = "DUP";
	public static final String DUP_TANDEM_STR = "DUP:TANDEM";
	public static final String DUP_INT_STR = "DUP:INT";
	public static final String INV_STR = "INV";
	public static final String INV_DUP_STR = "INVDUP";
	public static final String CNV_STR = "CNV";
	public static final String BND_STR = "BND";
	public static final String TRA_STR = "TRA";
    
	/**
     *
	 */
	public static final byte TYPE_INSERTION = 1;
	public static final byte TYPE_DELETION = 2;
	public static final byte TYPE_DEL_INV = 3;
	public static final byte TYPE_INVERSION = 4;
	public static final byte TYPE_INV_DUP = 5;
	public static final byte TYPE_DUPLICATION = 6;
	public static final byte TYPE_CNV = 7;
	public static final byte TYPE_BREAKEND = 8;
	public static final byte TYPE_TRANSLOCATION = 9;
    
	/**
	 * Confidence intervals of positions.
	 *
	 * Remark: some callers report a standard deviation instead of a confidence
	 * interval. Sniffles reports additional interval information in its BEDPE
	 * output (an alternative to VCF), which our programs disregard for 
	 * simplicity. Some callers use CILEN to express a "confidence interval 
	 * around inserted/deleted material between breakends": we interpret CILEN
	 * exactly like CIEND, and we ignore it for insertions (since representing
	 * variable-length insertion strings complicates our code).
	 */
	public static final String CI_SEPARATOR = ",";
	public static final String CIPOS_STR = "CIPOS";
	public static final String CIEND_STR = "CIEND";
	public static final String STD_START1_STR = "STD_quant_start";
	public static final String STD_START2_STR = "STD_POS1";
	public static final String STD_END1_STR = "STD_quant_stop";
	public static final String STD_END2_STR = "STD_POS2";
	public static final String CILEN_STR = "CILEN";
	public static final int CIPOS_STR_LENGTH = CIPOS_STR.length();
	public static final int CIEND_STR_LENGTH = CIEND_STR.length();
	public static final int STD_START1_STR_LENGTH = STD_START1_STR.length();
	public static final int STD_START2_STR_LENGTH = STD_START2_STR.length();
	public static final int STD_END1_STR_LENGTH = STD_END1_STR.length();
	public static final int STD_END2_STR_LENGTH = STD_END2_STR.length();
	public static final int CILEN_STR_LENGTH = CILEN_STR.length();
    
}