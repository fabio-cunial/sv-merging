import java.util.Arrays;
import java.io.*;


/**
 * Partitions a VCF into calls that overlap with a previously-computed list of 
 * intervals, and those that don't.
 */
public class FilterRegionsByIntervals {
    /**
     * INS are represented as intervals of constant length centered at the
     * insertion position.
     */
    private static final int INS_SLACK = 100;  // Arbitrary
    
    /**
     * Total lengths
     */
    private static final int CHR21_LENGTH = 46709983;
    private static final int CHR22_LENGTH = 50818468;
    
    /**
     * Rows: SV types (0=DEL, 1=INV, 2=DUP, 3=INS, 4=ALL).
     * Columns: positions on the chromosome.
     * Cell (i,j): number of calls of type $i$ that cover position $j$.
     */
    private static boolean[] keep_chr21, keep_chr22;
    
    
    /**
     * @param args 2: all and only the calls that have at least this fraction
     * of their bases marked in the regions file are assigned to the keep VCF.
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF = args[0];
        final String REGIONS_FILE = args[1];
        final double SV_FRACTION = Double.parseDouble(args[2]);
        final boolean ONLY_PASS = Integer.parseInt(args[3])==1;
        final String KEEP_VCF = args[4];
        final String DISCARD_VCF = args[5];
        
        loadRegions(REGIONS_FILE);
        filterVCF(INPUT_VCF,SV_FRACTION,ONLY_PASS,KEEP_VCF,DISCARD_VCF);
    }
    
    
    private static final void loadRegions(String inputFile) throws IOException {
        int i;
        int first, last;
        String str;
        BufferedReader br;
        boolean[] array;
        String[] tokens;
        
        keep_chr21 = new boolean[CHR21_LENGTH]; keep_chr22 = new boolean[CHR22_LENGTH];
        Arrays.fill(keep_chr21,false); Arrays.fill(keep_chr22,false);
        br = new BufferedReader(new FileReader(inputFile));
        str=br.readLine();
        while (str!=null) {
            tokens=str.split(",");
            first=Integer.parseInt(tokens[1]);
            last=Integer.parseInt(tokens[2]);
            if (tokens[0].equalsIgnoreCase("chr21")) array=keep_chr21;
            else if (tokens[0].equalsIgnoreCase("chr22")) array=keep_chr22;
            else array=null;
            if (array!=null) {
                for (i=first; i<=last; i++) array[i]=true;
            }
            str=br.readLine();
        }
        br.close();
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
     * Partitions $inputVcf$ in two files: the calls that have $>=svFraction$ 
     * fraction of their bases marked in a $keep_*$ array, and those that don't.
     */
    private static final void filterVCF(String inputVcf, double svFraction, boolean onlyPass, String outputVcfKeep, String outputVcfDiscard) throws IOException {
        int i;
        int length, chrLength, first, last, position, row, sum, to;
        String str;
        BufferedReader br;
        BufferedWriter bwKeep, bwDiscard;
        boolean[] keep;
        String[] tokens;
        int[][] histogram;
        
        bwKeep = new BufferedWriter(new FileWriter(outputVcfKeep));
        bwDiscard = new BufferedWriter(new FileWriter(outputVcfDiscard));
        br = new BufferedReader(new FileReader(inputVcf));
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)==COMMENT) {
                bwKeep.write(str); bwKeep.newLine();
                bwDiscard.write(str); bwDiscard.newLine();
                str=br.readLine();
                continue;
            }
            tokens=str.split("\t");
			if (onlyPass && !tokens[6].equalsIgnoreCase(PASS_STR)) {
				str=br.readLine();
				continue;
			}
            row=svType2Row(getField(tokens[7],SVTYPE_STR));
            if (row==-1) {
				str=br.readLine();
				continue;
            }
            position=Integer.parseInt(tokens[1]);
            if (row==3) { to=position+INS_SLACK; position-=INS_SLACK; }
            else {
                length=Integer.parseInt(getField(tokens[7],SVLEN_STR));
                if (length<0) length=-length;
                to=position+length-1;
            }
            if (tokens[0].equalsIgnoreCase("chr21")) {
                keep=keep_chr21;
                chrLength=CHR21_LENGTH;
            }
            else if (tokens[0].equalsIgnoreCase("chr22")) {
                keep=keep_chr22;
                chrLength=CHR22_LENGTH;
            }
            else { keep=null; chrLength=0; }
            if (to>=chrLength) to=chrLength-1;
            sum=0;
            for (i=position; i<=to; i++) {
                if (keep[i]) sum++;
            }
            if (sum>=(to-position+1)*svFraction) { bwKeep.write(str); bwKeep.newLine(); }
            else { bwDiscard.write(str); bwDiscard.newLine(); }
            str=br.readLine();
        }
        br.close(); bwKeep.close(); bwDiscard.close();
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