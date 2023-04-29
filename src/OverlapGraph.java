import java.util.Arrays;
import java.io.*;


/**
 * Prints basic statistics on the graph where nodes are calls and edges are
 * overlaps.
 */
public class OverlapGraph {
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
    private static int[] componentSize;
    private static double[] componentTrfOverlap;  // Fraction, avg. on all calls
    private static double[][] endpointAvg, endpointVariance;
    private static double[] lengthAvg, lengthVariance;
    private static int[][] endpointMin, endpointMax;
    private static int[] lengthMin, lengthMax;
    private static int lastComponent;
    
    /**
     * Columns: positions on the chromosome.
     */
    private static boolean[] trf_mask_chr21, trf_mask_chr22;
    
    /**
     * Output histograms
     */
    private static int[] histogram_degree, histogram_componentSize;
    
    
    /**
     * @param args 3: list of zero-based column IDs of the VCF matrix 
     * (=individuals); set it to $null$ to plot overlaps of all SVs in
     * $VCF_FILE$.
     */
    public static void main(String[] args) throws Exception {
        final String VCF_FILE = args[0];
        final String ONLY_CHR = args[1];
        final int SV_TYPE = Integer.parseInt(args[2]);
        final String SAMPLE_ID_LIST = args[3];
        final boolean ONLY_PASS = Integer.parseInt(args[4])==1;
        final int MAX_DEGREE = Integer.parseInt(args[5]);
        final int MAX_COMPONENT_SIZE = Integer.parseInt(args[6]);
        final String TRF_FILE = args[7];
        final String OUTPUT_DIR = args[8];
        
        int sampleID;
        String str;
        BufferedReader br;
        
        trf_mask_chr21 = new boolean[CHR21_LENGTH];
        trf_mask_chr22 = new boolean[CHR22_LENGTH];
        buildTrfMasks(TRF_FILE);
        if (SAMPLE_ID_LIST.equalsIgnoreCase("null")) {
            System.err.print("Loading all calls... ");
            loadCalls(-1,ONLY_CHR,SV_TYPE,VCF_FILE,ONLY_PASS);
            System.err.println("DONE. "+(lastCall+1)+" calls loaded.");
            System.err.print("Building connected components... ");
            buildComponents(MAX_DEGREE,MAX_COMPONENT_SIZE);
            System.err.println("DONE. "+(lastComponent+1)+" components.");
            printHistogram(-1,histogram_degree,ONLY_CHR+"_"+SV_TYPES[SV_TYPE]+"_degree_histogram",OUTPUT_DIR);
            printHistogram(-1,histogram_componentSize,ONLY_CHR+"_"+SV_TYPES[SV_TYPE]+"_componentSize_histogram",OUTPUT_DIR);
            printGraph(OUTPUT_DIR+"/all_"+ONLY_CHR+"_"+SV_TYPES[SV_TYPE]+"_graph.dot");
            printComponentStats(OUTPUT_DIR+"/all_"+ONLY_CHR+"_"+SV_TYPES[SV_TYPE]+"_componentStats.txt");
            printCallTrfOverlap(OUTPUT_DIR+"/all_"+ONLY_CHR+"_"+SV_TYPES[SV_TYPE]+"_trfOverlaps.txt");
        }
        else {
            br = new BufferedReader(new FileReader(SAMPLE_ID_LIST));
            str=br.readLine();
            while (str!=null) {
                sampleID=Integer.parseInt(str);
                System.err.print("Loading calls in sample "+sampleID+"... ");
                loadCalls(sampleID,ONLY_CHR,SV_TYPE,VCF_FILE,ONLY_PASS);
                System.err.println("DONE. "+(lastCall+1)+" calls loaded.");
                System.err.print("Building connected components... ");
                buildComponents(MAX_DEGREE,MAX_COMPONENT_SIZE);
                System.err.println("DONE. "+(lastComponent+1)+" components.");
                printHistogram(sampleID,histogram_degree,ONLY_CHR+"_"+SV_TYPES[SV_TYPE]+"_degree_histogram",OUTPUT_DIR);
                printHistogram(sampleID,histogram_componentSize,ONLY_CHR+"_"+SV_TYPES[SV_TYPE]+"_componentSize_histogram",OUTPUT_DIR);
                str=br.readLine();
            }
            br.close();
        }
    }
    
    
	private static final void buildTrfMasks(String path) throws IOException {
		int i;
		int start, end;
		String str;
		BufferedReader br;
        boolean[] mask;
		String[] tokens;

		br = new BufferedReader(new FileReader(path));
		str=br.readLine();
		while (str!=null) {
			tokens=str.split(",");
            start=Integer.parseInt(tokens[1])-1;
            if (start<0) start=0;
            end=Integer.parseInt(tokens[2])-1;
            mask=null;
            if (tokens[0].equalsIgnoreCase("chr21")) {
                if (end>=CHR21_LENGTH) end=CHR21_LENGTH-1;
                mask=trf_mask_chr21;
            }
            else if (tokens[0].equalsIgnoreCase("chr22")) {
                if (end>=CHR22_LENGTH) end=CHR22_LENGTH-1;
                mask=trf_mask_chr22;
            }
            if (mask!=null) {
                for (i=start; i<=end; i++) mask[i]=true;
            }
			str=br.readLine();
		}
		br.close();
    }
    
    
    /**
     * Remark: the elements of $calls$ are not necessarily sorted when the
     * procedure completes.
     *
     * @param sampleID zero-based ID of a column of the VCF matrix 
     * (=individual); -1=build the graph of all calls, without checking if a
     * call is active in a specific individual;
     * @param svType load only calls of a specific type: 0=DEL, 1=INV, 2=DUP,
     * 3=INS, 4=ALL.
     */
    private static final void loadCalls(int sampleID, String onlyChr, int svType, String path, boolean onlyPass) throws IOException {
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
            if (!tokens[0].equalsIgnoreCase(onlyChr)) {
				str=br.readLine();
				continue;
            }
            position=Integer.parseInt(tokens[1]);
			if (onlyPass && !tokens[6].equalsIgnoreCase(PASS_STR)) {
				str=br.readLine();
				continue;
			}
            row=svType2Row(getField(tokens[7],SVTYPE_STR));
            if ((svType==4 && row==-1) || (svType!=4 && row!=svType)) {
				str=br.readLine();
				continue;
            }
            if (sampleID>=0) n=(tokens[9+sampleID].charAt(0)=='1'?1:0)+(tokens[9+sampleID].charAt(2)=='1'?1:0);
            else n=1;
            if (n!=0) {
                if (row==3) { to=position+INS_SLACK; position-=INS_SLACK; }
                else {
                    length=Integer.parseInt(getField(tokens[7],SVLEN_STR));
                    if (length<0) length=-length;
                    to=position+length-1;
                }
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
            }
            str=br.readLine();
        }
        br.close();
    }
    
    
    /**
     * Computes degree and connected component of every call in $calls$, and 
     * builds $histogram_degree, histogram_componentSize$.
     */
    private static final void buildComponents(int maxDegree, int maxComponentSize) {
        final int N_BINS = 50;  // Arbitrary
        int i, j;
        int degree, size, chr, first, last, component, length;
        double quantum;
        
        System.err.print("Sorting calls...");
        Arrays.sort(calls,0,lastCall+1);
        System.err.println("DONE");
        histogram_degree = new int[maxDegree+1];
        histogram_componentSize = new int[maxComponentSize+1];
        lastComponent=-1;
        for (i=0; i<=lastCall; i++) {
            chr=calls[i].chr; first=calls[i].first; last=calls[i].last;
            if (calls[i].component==-1) calls[i].component=++lastComponent;
            for (j=i+1; j<=lastCall; j++) {
                if (calls[j].chr!=chr || calls[j].first>last) break;
                calls[j].component=calls[i].component;
                calls[i].degree++; calls[j].degree++;
            }
        }
        componentSize = new int[lastComponent+1];
        componentTrfOverlap = new double[lastComponent+1];
        endpointAvg = new double[lastComponent+1][2];
        lengthAvg = new double[lastComponent+1];
        endpointMin = new int[lastComponent+1][2];
        for (i=0; i<=lastComponent; i++) { endpointMin[i][0]=Integer.MAX_VALUE; endpointMin[i][1]=Integer.MAX_VALUE; }
        endpointMax = new int[lastComponent+1][2];
        lengthMin = new int[lastComponent+1];
        for (i=0; i<=lastComponent; i++) lengthMin[i]=Integer.MAX_VALUE;
        lengthMax = new int[lastComponent+1];
        for (i=0; i<=lastCall; i++) {
            degree=calls[i].degree;
            if (degree>maxDegree) degree=maxDegree;
            histogram_degree[degree]++;
            component=calls[i].component;
            componentSize[component]++;
            componentTrfOverlap[component]+=calls[i].getTrfOverlap();
            if (calls[i].first<endpointMin[component][0]) endpointMin[component][0]=calls[i].first;
            if (calls[i].first>endpointMax[component][0]) endpointMax[component][0]=calls[i].first;
            endpointAvg[component][0]+=calls[i].first;
            if (calls[i].last<endpointMin[component][1]) endpointMin[component][1]=calls[i].last;
            if (calls[i].last>endpointMax[component][1]) endpointMax[component][1]=calls[i].last;
            endpointAvg[component][1]+=calls[i].last;
            length=calls[i].last-calls[i].first+1;
            if (length<lengthMin[component]) lengthMin[component]=length;
            if (length>lengthMax[component]) lengthMax[component]=length;
            lengthAvg[component]+=length;
        }
        for (i=0; i<=lastComponent; i++) {
            size=componentSize[i];
            componentTrfOverlap[i]/=size;
            endpointAvg[i][0]/=size;
            endpointAvg[i][1]/=size;
            lengthAvg[i]/=size;
            if (size>maxComponentSize) size=maxComponentSize;
            histogram_componentSize[size]++;
        }
        endpointVariance = new double[lastComponent+1][2];
        lengthVariance = new double[lastComponent+1];
        for (i=0; i<=lastCall; i++) {
            component=calls[i].component;
            endpointVariance[component][0]+=(calls[i].first-endpointAvg[component][0])*(calls[i].first-endpointAvg[component][0]);
            endpointVariance[component][1]+=(calls[i].last-endpointAvg[component][1])*(calls[i].last-endpointAvg[component][1]);
            length=calls[i].last-calls[i].first+1;
            lengthVariance[component]+=(length-lengthAvg[component])*(length-lengthAvg[component]);
        }
        for (i=0; i<=lastComponent; i++) {
            size=componentSize[i];
            endpointVariance[i][0]/=size;
            endpointVariance[i][1]/=size;
            lengthVariance[i]/=size;
        }
    }
    
    
    /**
     *
     */
    private static final void printComponentStats(String path) throws IOException {
        int i;
        BufferedWriter bw;
        
        bw = new BufferedWriter(new FileWriter(path));
        for (i=0; i<=lastComponent; i++) bw.write(componentSize[i]+","+componentTrfOverlap[i]+","+endpointVariance[i][0]+","+endpointVariance[i][1]+","+lengthVariance[i]+","+(endpointMax[i][0]-endpointMin[i][0])+","+(endpointMax[i][1]-endpointMin[i][1])+","+(lengthMax[i]-lengthMin[i])+"\n");
        bw.close();
    }
    
    
    /**
     * Remark: the procedure assumes that $calls$ has already been sorted.
     */
    private static final void printGraph(String path) throws IOException {
        int i, j;
        int chr, first, last;
        BufferedWriter bw;
        
        bw = new BufferedWriter(new FileWriter(path));
        bw.write("graph G {");
        for (i=0; i<=lastCall; i++) {
            chr=calls[i].chr; first=calls[i].first; last=calls[i].last;
            bw.write(i+" "+calls[i].toDot()+";\n");
            for (j=i+1; j<=lastCall; j++) {
                if (calls[j].chr!=chr || calls[j].first>last) break;
                bw.write(i+" -- "+j+";\n");
            }
        }
        bw.write("}"); bw.close();
    }
    
    
    /**
     * Prints $histogram$ to a file, prepended by $sampleID$.
     *
     * @param sampleID zero-based; -1=all individuals.
     */
    private static final void printHistogram(int sampleID, int[] histogram, String histogramName, String outputDir) throws IOException {
        int i;
        BufferedWriter bw;
        
        bw = new BufferedWriter(new FileWriter(outputDir+"/"+(sampleID>=0?sampleID:"all")+"_"+histogramName+".txt"));
        bw.write(sampleID+"");
        for (i=0; i<histogram.length; i++) bw.write(","+histogram[i]);
        bw.newLine();
        bw.close();
    }
    
    
    /**
     * Prints a row for every call.
     */
    private static final void printCallTrfOverlap(String path) throws IOException {
        int i;
        BufferedWriter bw;
        
        bw = new BufferedWriter(new FileWriter(path));
        for (i=0; i<=lastCall; i++) bw.write((calls[i].last-calls[i].first+1)+","+calls[i].getTrfOverlap()+"\n");
        bw.close();
    }
    
    
    private static class Call implements Comparable {
        public int type, chr, first, last;
        public int component, degree;
        
        public Call() { }
        
        public Call(int t, int c, int f, int l) {
            this.type=t; this.chr=c; this.first=f; this.last=l;
            component=-1; degree=0;
        }
        
        /**
         * @return the fraction of the call that overlaps with the TRF mask.
         */
        public double getTrfOverlap() {
            int i;
            double out;
            final boolean[] mask = chr==21?trf_mask_chr21:trf_mask_chr22;
            
            out=0.0;
            for (i=first; i<=last; i++) {
                if (mask[i]) out+=1.0;
            }
            return out/(last-first+1);
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
            return chr+"\t"+first+"\t"+SV_TYPES[type]+"\t"+last; 
        }
        
        public String toDot() {
            return "[type=\""+SV_TYPES[type]+"\",chr=\""+chr+"\",first=\""+first+"\",last=\""+last+"\",degree=\""+degree+"\",component=\""+component+"\"]";
        }
    }
    
    
	private static final int svType2Row(String type) {
		if (type==null || type.length()==0) return -1;
		if ( type.equalsIgnoreCase(DEL_STR) || 
			 type.equalsIgnoreCase(DEL_ME_STR)
		   ) return 0;
		else if (type.equalsIgnoreCase(INV_STR)) return 1;
        else if ( type.equalsIgnoreCase(DUP_STR) ||
			      type.equalsIgnoreCase(DUP_TANDEM_STR) ||
				  type.equalsIgnoreCase(DUP_INT_STR)
			    ) return 2;
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