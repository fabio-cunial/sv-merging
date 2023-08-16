import java.util.Arrays;
import java.io.*;


/**
 * On jasmine, ~50% of all calls <=10kb are discarded because they overlap
 * another call. ~80% of all calls >10kb are discarded because they overlap
 * another call.
 */
public class KeepSimpleCalls {
    /**
     * Calls and tandems
     */
    private static Call[] calls, calls_short, calls_long;
    private static int lastCall, lastCall_short, lastCall_short_keep, lastCall_long, lastCall_long_keep;
    
    /**
     * VCF header
     */
    private static StringBuilder header;
    
    
    /**
     * 
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF = args[0];
        final boolean PASS_ONLY = Integer.parseInt(args[1])==1;
        final String TRF_FILE = args[2];
        final int IDENTITY_THRESHOLD = Integer.parseInt(args[3]);
        final int SV_LENGTH_THRESHOLD = Integer.parseInt(args[4]);
        final String OUTPUT_VCF_SHORT = args[5];
        final String OUTPUT_VCF_LONG = args[6];
        
        int i;
        BufferedReader br;
        BufferedWriter bw;
        
        // Loading calls and tandems
        System.err.println("Loading all calls...");
        loadCalls(INPUT_VCF,PASS_ONLY);
        System.err.println("DONE. "+(lastCall+1)+" calls loaded.");
        splitByLength(SV_LENGTH_THRESHOLD);
        System.err.println("Loading tandems...");
        loadTandems(TRF_FILE);
        System.err.println("DONE. "+(lastCall+1)+" tandems loaded.");
        duplicateTandems();
        
        // Filtering calls
        System.err.println("Filtering calls... ");
        lastCall_short_keep=filterCalls(true,IDENTITY_THRESHOLD);
        lastCall_long_keep=filterCalls(false,IDENTITY_THRESHOLD);
        
        // Outputting
        bw = new BufferedWriter(new FileWriter(OUTPUT_VCF_SHORT));
        bw.write(header.toString());
        for (i=0; i<=lastCall_short_keep; i++) bw.write(calls_short[i].str+"\n");
        bw.close();
        bw = new BufferedWriter(new FileWriter(OUTPUT_VCF_LONG));
        bw.write(header.toString());
        for (i=0; i<=lastCall_long_keep; i++) bw.write(calls_long[i].str+"\n");
        bw.close();
    }
    
    
    /**
     * Loads the calls in $path$ into global array $calls$, which is not
     * necessarily sorted when the procedure completes.
     */
    private static final void loadCalls(String path, boolean passOnly) throws IOException {
        final int INS_SLACK = 100;  // Arbitrary
        final int CAPACITY = 1000;  // Arbitrary
        int i;
        int position, length, row, to, chr;
        String str;
        BufferedReader br;
        String[] tokens;
        
        calls = new Call[CAPACITY]; lastCall=-1;
        header = new StringBuilder();
        br = new BufferedReader(new FileReader(path));
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)==COMMENT) {
                header.append(str+"\n");
                str=br.readLine();
                continue;
            }
            tokens=str.split("\t");
            position=Integer.parseInt(tokens[1]);
			if (passOnly && !tokens[6].equalsIgnoreCase(PASS_STR) && !tokens[6].equalsIgnoreCase(".")) {
				str=br.readLine();
				continue;
			}
            row=svType2Row(getField(tokens[7],SVTYPE_STR));
            if (row==-1) {
				str=br.readLine();
				continue;
            }
            if (row==3) { to=position+INS_SLACK; position-=INS_SLACK; }
            else {
                length=Integer.parseInt(getField(tokens[7],SVLEN_STR));
                if (length<0) length=-length;
                to=position+length;
            }
            chr=chr2int(tokens[0]);
            if (chr==-1) {
                str=br.readLine();
                continue;
            }
            lastCall++;
            if (lastCall==calls.length) {
                Call[] newArray = new Call[calls.length<<1];
                System.arraycopy(calls,0,newArray,0,calls.length);
                calls=newArray;
            }
            calls[lastCall] = new Call(row,chr,position+1,to,str);
            str=br.readLine();
        }
        br.close();
    }
    
    
    /**
     * Partitions the calls in $calls$ by length into $calls_short$ and 
     * $calls_long$, using $threshold$. At the end of the procedure, all
     * elements in $calls$ are NULL, but $calls$ itself is not null.
     */
    private static final void splitByLength(int threshold) {
        int i, j;
        Call tmpCall;
        
        lastCall_short=-1;
        for (i=0; i<=lastCall; i++) {
            if (calls[i].last-calls[i].first+1<threshold) {
                lastCall_short++;
                tmpCall=calls[lastCall_short];
                calls[lastCall_short]=calls[i];
                calls[i]=tmpCall;
            }
        }
        calls_short = new Call[lastCall_short+1];
        System.arraycopy(calls,0,calls_short,0,lastCall_short+1);
        lastCall_long=lastCall-lastCall_short-1;
        calls_long = new Call[lastCall_long+1];
        System.arraycopy(calls,lastCall_short+1,calls_long,0,lastCall_long+1);
        for (i=0; i<=lastCall; i++) calls[i]=null;
        lastCall=-1;
    }
    
    
    /**
     * Tandem intervals are modeled as calls of a special type, and they are
     * loaded in $calls$ as well.
     */
	private static final void loadTandems(String path) throws IOException {
		int i;
		int chr, start, end;
		String str;
		BufferedReader br;
		String[] tokens;
        
		br = new BufferedReader(new FileReader(path));
		str=br.readLine();
		while (str!=null) {
			tokens=str.split("\t");
            chr=chr2int(tokens[0]);
            if (chr==-1) {
                str=br.readLine();
                continue;
            }
            start=Integer.parseInt(tokens[1]);
            end=Integer.parseInt(tokens[2]);
            lastCall++;
            if (lastCall==calls.length) {
                Call[] newArray = new Call[calls.length<<1];
                System.arraycopy(calls,0,newArray,0,calls.length);
                calls=newArray;
            }
            calls[lastCall] = new Call(-1,chr,start,end,null);
			str=br.readLine();
		}
		br.close();
    }
    
    
    /**
     * Copies all the elements in $calls$ to the end of $calls_short$ and 
     * $calls_long$. At the end of the procedure $calls$ is NULL.
     */
    private static final void duplicateTandems() {
        int i;
        
        if (lastCall_short+1+lastCall+1>=calls_short.length) {
            Call[] newArray = new Call[lastCall_short+1+lastCall+1];
            System.arraycopy(calls_short,0,newArray,0,lastCall_short+1);
            System.arraycopy(calls,0,newArray,lastCall_short+1,lastCall+1);
        }
        if (lastCall_long+1+lastCall+1>=calls_long.length) {
            Call[] newArray = new Call[lastCall_long+1+lastCall+1];
            System.arraycopy(calls_long,0,newArray,0,lastCall_long+1);
            System.arraycopy(calls,0,newArray,lastCall_long+1,lastCall+1);
        }
        for (i=0; i<=lastCall; i++) calls[i]=null;
        calls=null; lastCall=-1;
    }
    
    
    /**
     * Discards calls that are not isolated or that overlap tandem repeats.
     * 
     * @param shortOrLong TRUE=short, FALSE=long.
     * @return X, where array [0..X] contains all the calls that were kept and
     * [X+1..] contains all the calls that were discarded at the end of the
     * procedure.
     */
    private static final int filterCalls(boolean shortOrLong, int identityThreshold) {
        final Call[] calls = shortOrLong?calls_short:calls_long;
        final int lastCall = shortOrLong?lastCall_short:lastCall_long;
        int i, j;
        int chr, first, last, lastComponent;
        int nFilteredOut_componentSize, nFilteredOut_tandem;
        Call tmpCall;
        int[] componentSize_calls;
        
        System.err.println("Sorting calls...");
        Arrays.sort(calls,0,lastCall+1);
        System.err.println("DONE");
        System.err.println("Computing connected components...");
        lastComponent=-1;
        for (i=0; i<=lastCall; i++) {
            chr=calls[i].chr; first=calls[i].first; last=calls[i].last;
            if (calls[i].component==-1) calls[i].component=++lastComponent;
            for (j=i+1; j<=lastCall; j++) {
                if (calls[j].chr!=chr || calls[j].first>last+identityThreshold) break;
                calls[j].component=calls[i].component;
            }
        }
        System.err.println("Filtering calls...");
        componentSize_calls = new int[lastComponent+1];
        for (i=0; i<=lastCall; i++) {
            if (calls[i].type!=-1) componentSize_calls[calls[i].component]++;
        }
        j=-1; nFilteredOut_componentSize=0; nFilteredOut_tandem=0;
        for (i=0; i<=lastCall; i++) {
            if (calls[i].type==-1) continue;
            if (componentSize_calls[calls[i].component]>1) { nFilteredOut_componentSize++; continue; }
            if (inTandem(i,calls,lastCall,identityThreshold)) { nFilteredOut_tandem++; continue; }
            j++; tmpCall=calls[j]; calls[j]=calls[i]; calls[i]=tmpCall;
        }
        System.err.println("Discarded "+nFilteredOut_componentSize+" "+(shortOrLong?"short":"long")+" calls because they are not isolated ("+((100.0*nFilteredOut_componentSize)/(lastCall+1))+"%).");
        System.err.println("Discarded "+nFilteredOut_tandem+" "+(shortOrLong?"short":"long")+" calls because they straddle or are contained in a tandem.");
        return j;
    }
    
    
    /**
     * @return TRUE iff $calls[callID]$ is contained in or straddles a tandem.
     */
    private static final boolean inTandem(int callID, Call[] calls, int lastCall, int identityThreshold) {
        int i;
        final Call call = calls[callID];
        final int chr = call.chr;
        final int first = call.first;
        final int last = call.last;
        final int component = call.component;
        
        for (i=callID+1; i<=lastCall; i++) {
            if (calls[i].chr!=chr || calls[i].first>last+identityThreshold) break;
            if (calls[i].component!=component || calls[i].type!=-1) continue;
            if ( (first>=calls[i].first-identityThreshold && last<=calls[i].last+identityThreshold) ||
                 (last>=calls[i].first && last<=calls[i].last)
               ) return true;
        }
        for (i=callID-1; i>=0; i--) {
            if (calls[i].chr!=chr) break;
            if (calls[i].component!=component || calls[i].type!=-1) continue;
            if ( (first>=calls[i].first-identityThreshold && last<=calls[i].last+identityThreshold) ||
                 (first>=calls[i].first && first<=calls[i].last)
               ) return true;
        }
        return false;
    }
    
    
    /**
     * @return the one-based integer associated with the chromosome.
     */
    private static final int chr2int(String chr) {
        if (chr.length()>5) return -1;
        final char c = chr.charAt(3);
        if (Character.isDigit(c)) return Integer.parseInt(chr.substring(3));
        switch (c) {
            case 'x': return 23;
            case 'X': return 23;
            case 'y': return 24;
            case 'Y': return 24;
            case 'm': return 25;
            case 'M': return 25;
            default: return -1;
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
    
    
    private static class Call implements Comparable {
        public int chr;  // One-based
        public int first, last;  // One-based
        public int type;  // -1=tandem interval.
        public String str;
        
        public int component;
        
        public Call() { }
        
        public Call(int t, int c, int f, int l, String s) {
            this.type=t; this.chr=c; this.first=f; this.last=l;
            this.str=s; component=-1;
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