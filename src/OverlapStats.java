import java.util.Arrays;
import java.io.*;


public class OverlapStats {
    /**
     * Total lengths
     */
    private static final int CHR21_LENGTH = 46709983;
    private static final int CHR22_LENGTH = 50818468;
    
    /**
     * Rows: SV types (0=DEL, 1=INV, 2=DUP, 3=INS).
     * Columns: positions on the chromosome.
     */
    private static int[][] histogram_chr21, histogram_chr22;
    private static int n00, n01, n11, nDD, nD0, nD1;
    
    
    public static void main(String[] args) throws Exception {
        final String VCF_FILE = args[0];
        final String SAMPLE_ID_LIST = args[1];
        final boolean ONLY_PASS = Integer.parseInt(args[2])==1;
        final String OUTPUT_DIR = args[3];
        final int MAX_OVERLAP = Integer.parseInt(args[4]);
        
        int sampleID;
        String str;
        BufferedReader br;
        BufferedWriter bw;
        
        histogram_chr21 = new int[4][CHR21_LENGTH];
        histogram_chr22 = new int[4][CHR22_LENGTH];
        br = new BufferedReader(new FileReader(SAMPLE_ID_LIST));
        str=br.readLine();
        while (str!=null) {
            sampleID=Integer.parseInt(str);
            buildHistograms(VCF_FILE,sampleID,ONLY_PASS);
            printHistogram(OUTPUT_DIR,21,0,sampleID,MAX_OVERLAP);
            printHistogram(OUTPUT_DIR,21,1,sampleID,MAX_OVERLAP);
            printHistogram(OUTPUT_DIR,21,2,sampleID,MAX_OVERLAP);
            printHistogram(OUTPUT_DIR,21,3,sampleID,MAX_OVERLAP);
            printHistogram(OUTPUT_DIR,22,0,sampleID,MAX_OVERLAP);
            printHistogram(OUTPUT_DIR,22,1,sampleID,MAX_OVERLAP);
            printHistogram(OUTPUT_DIR,22,2,sampleID,MAX_OVERLAP);
            printHistogram(OUTPUT_DIR,22,3,sampleID,MAX_OVERLAP);
            bw = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/"+sampleID+"_gtCounts.txt"));
            bw.write(sampleID+","+n00+","+n01+","+n11+","+nDD+","+nD0+","+nD1+"\n");
            bw.close();
            str=br.readLine();
        }
        br.close();
    }
    
    
    /**
     * Updates variables $histogram_*$ and $n*$.
     */
    private static final void buildHistograms(String path, int sampleID, boolean onlyPass) throws IOException {
        int i, n;
        int position, length, row, to;
        String str, gt;
        BufferedReader br;
        String[] tokens;
        
        n00=0; n01=0; n11=0; nDD=0; nD0=0; nD1=0;
        for (i=0; i<histogram_chr21.length; i++) Arrays.fill(histogram_chr21[i],0);
        for (i=0; i<histogram_chr22.length; i++) Arrays.fill(histogram_chr22[i],0);
        br = new BufferedReader(new FileReader(path));
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)==COMMENT) {
                str=br.readLine();
                continue;
            }
            tokens=str.split("\t");
            position=Integer.parseInt(tokens[1]);
			if (onlyPass && !tokens[6].equalsIgnoreCase(PASS_STR)) {
				str=br.readLine();
				continue;
			}
            gt=tokens[9+sampleID].substring(0,3);
            if (gt.equals("0/0")) n00++;
            else if (gt.equals("0/1") || gt.equals("1/0")) n01++;
            else if (gt.equals("1/1")) n11++;
            else if (gt.equals("./.")) nDD++;
            else if (gt.equals("./0") || gt.equals("0/.")) nD0++;
            else if (gt.equals("./1") || gt.equals("1/.")) nD1++;
            row=svType2Row(getField(tokens[7],SVTYPE_STR));     
            if (row==-1) {
				str=br.readLine();
				continue;
            }
            length=Integer.parseInt(getField(tokens[7],SVLEN_STR));
            if (length<0) length=-length;
            n=(tokens[9+sampleID].charAt(0)=='1'?1:0)+(tokens[9+sampleID].charAt(2)=='1'?1:0);
            if (n!=0) {
                to=position+length-1;
                if (tokens[0].equals("chr21")) {
                    if (to>=CHR21_LENGTH) to=CHR21_LENGTH-1;
                    for (i=position-1; i<=to; i++) histogram_chr21[row][i]++;
                }
                else if (tokens[0].equals("chr22")) {
                    if (to>=CHR22_LENGTH) to=CHR22_LENGTH-1;
                    for (i=position-1; i<=to; i++) histogram_chr22[row][i]++;
                }
            }
            str=br.readLine();
        }
        br.close();
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
     * Prints a file with a single row that contains the probability that a 
     * position of the genome covered by an SV, is covered by X total SVs, for
     * $X \in [1..maxOverlap]$.
     */
    private static final void printHistogram(String outputDir, int chr, int row, int sampleID, int maxOverlap) throws IOException {
        int i, n;
        int length;
        long sum;
        String svType;
        BufferedWriter bw;
        int[][] histogram;
        double[] overlapHistogram;
        
        overlapHistogram = new double[maxOverlap+1];
        if (chr==21) { histogram=histogram_chr21; length=CHR21_LENGTH; }
        else if (chr==22) { histogram=histogram_chr22; length=CHR22_LENGTH; }
        else { histogram=null; length=-1; }
        sum=0;
        for (i=0; i<length; i++) {
            if (histogram[row][i]==0) continue;
            sum++;
            n=histogram[row][i]>maxOverlap?maxOverlap:histogram[row][i];
            overlapHistogram[n]++;
        }
        svType="";
        switch (row) {
            case 0: svType="del"; break;
            case 1: svType="inv"; break;
            case 2: svType="dup"; break;
            case 3: svType="ins"; break;
        }
        bw = new BufferedWriter(new FileWriter(outputDir+"/"+sampleID+"_chr"+chr+"_"+svType+"_histogram.txt"));
        bw.write(sampleID+",");
        for (i=1; i<=maxOverlap; i++) bw.write((overlapHistogram[i]/sum)+",");
        bw.newLine(); bw.close();
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