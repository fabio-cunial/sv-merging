import java.util.Arrays;
import java.io.*;


/**
 * Prints in output all maximal low-overlap regions of the chromosomes.
 */
public class KeepLowOverlapRegions {
    /**
     * SV types. INS are represented as intervals of constant length centered at
     * the insertion position.
     */
    private static final String[] SV_TYPES = new String[] {"del","inv","dup","ins","all"};
    private static final int N_SV_TYPES = 5;
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
    private static int[][] histogram_chr21, histogram_chr22;
    
    /**
     * Columns: positions on the chromosome.
     */
    private static boolean[] trf_mask_chr21, trf_mask_chr22;
    
    
    /**
     * 
     */
    public static void main(String[] args) throws IOException {
        final String VCF_FILE = args[0];
        final int MAX_OVERLAP = Integer.parseInt(args[1]);
        final boolean ONLY_PASS = Integer.parseInt(args[2])==1;
        final String TRF_FILE = args[3];
        final String OUTPUT_DIR = args[4];
        
        BufferedWriter bw;
        
        // Building overlap counts per position
        histogram_chr21 = new int[5][CHR21_LENGTH];
        histogram_chr22 = new int[5][CHR22_LENGTH];
        buildPositionOverlapHistograms(VCF_FILE,ONLY_PASS);
        
        // Building TRF masks
        trf_mask_chr21 = new boolean[CHR21_LENGTH];
        trf_mask_chr22 = new boolean[CHR22_LENGTH];
        buildTrfMasks(TRF_FILE);
        
        // Computing basic statistics on the intersection high-overlap / TRF
        highOverlapIntersectTrf(21,MAX_OVERLAP+1);
        highOverlapIntersectTrf(22,MAX_OVERLAP+1);
        
        // Printing maximal low-overlap regions
        bw = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/keep_"+MAX_OVERLAP+".txt"));
        keepRegions(21,N_SV_TYPES-1,MAX_OVERLAP,bw);
        keepRegions(22,N_SV_TYPES-1,MAX_OVERLAP,bw);
        bw.close();
        
        // Statistics: regions with maximum overlap of any type.
        bw = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/maxOverlap.txt"));
        maxOverlapRegions(21,N_SV_TYPES-1,2,bw);
        maxOverlapRegions(22,N_SV_TYPES-1,2,bw);
        bw.close();
        
        // Statistics: regions with maximum overlap of any type, outside TRF.
        bw = new BufferedWriter(new FileWriter(OUTPUT_DIR+"/maxOverlapNoTRF.txt"));
        maxOverlapRegions(21,N_SV_TYPES-1,0,bw);
        maxOverlapRegions(22,N_SV_TYPES-1,0,bw);
        bw.close();
    }
    
    
    /**
     * Updates variables $histogram_*$ using all SVs, without focusing on a
     * specific individual.
     */
    private static final void buildPositionOverlapHistograms(String path, boolean onlyPass) throws IOException {
        int i;
        int position, length, row, to;
        String str;
        BufferedReader br;
        String[] tokens;
        
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
            if (tokens[0].equals("chr21")) {
                if (to>=CHR21_LENGTH) to=CHR21_LENGTH-1;
                for (i=position; i<=to; i++) histogram_chr21[row][i]++;
                for (i=position; i<=to; i++) histogram_chr21[N_SV_TYPES-1][i]++;
            }
            else if (tokens[0].equals("chr22")) {
                if (to>=CHR22_LENGTH) to=CHR22_LENGTH-1;
                for (i=position; i<=to; i++) histogram_chr22[row][i]++;
                for (i=position; i<=to; i++) histogram_chr22[N_SV_TYPES-1][i]++;
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
     * Appends to $bw$ every maximal interval of $chr$ that is covered by 
     * $<=maxValue$ calls of type $svType$.
     *
     * @param svType a row of a $histogram_*$ matrix.
     */
    private static final void keepRegions(int chr, int svType, int maxValue, BufferedWriter bw) throws IOException {
        int i;
        int chrLength, first, last;
        int[][] histogram;
        
        if (chr==21) { histogram=histogram_chr21; chrLength=CHR21_LENGTH; }
        else { histogram=histogram_chr22; chrLength=CHR22_LENGTH; }
        first=-1; last=-1;
        for (i=0; i<chrLength; i++) {
            if (histogram[svType][i]<=maxValue) {
                if (first==-1) first=i;
                last=i;
            }
            else {
                if (first!=-1) bw.write("chr"+chr+","+first+","+last+"\n");
                first=-1; last=-1;
            }
        }
        if (first!=-1) bw.write("chr"+chr+","+first+","+last+"\n");
    }
    
    
    /**
     * Appends to $bw$ every maximal interval of $chr$ that is covered by the
     * maximum possible number of calls of type $svType$.
     *
     * @param svType a row of a $histogram_*$ matrix;
     * @param trfMode 0=do not consider positions inside the TRF mask; 
     * 1=consider only positions inside the TRF mask; 2=do not consider the TRF
     * mask at all.
     */
    private static final void maxOverlapRegions(int chr, int svType, int trfMode, BufferedWriter bw) throws IOException {
        int i;
        int chrLength, first, last, maxValue;
        boolean[] mask;
        int[][] histogram;
        
        if (chr==21) { histogram=histogram_chr21; chrLength=CHR21_LENGTH; mask=trf_mask_chr21; }
        else { histogram=histogram_chr22; chrLength=CHR22_LENGTH; mask=trf_mask_chr22; }
        maxValue=0;
        for (i=0; i<chrLength; i++) {
            if ((trfMode==0 && mask[i]) || (trfMode==1 && !mask[i])) continue;
            if (histogram[svType][i]>maxValue) maxValue=histogram[svType][i];
        }
        System.err.println("Max overlap for chr"+chr+": "+maxValue);
        first=-1; last=-1;
        for (i=0; i<chrLength; i++) {
            if ( ((trfMode==0 && !mask[i]) || (trfMode==1 && mask[i]) || trfMode==2) && histogram[svType][i]==maxValue ) {
                if (first==-1) first=i;
                last=i;
            }
            else {
                if (first!=-1) bw.write("chr"+chr+","+first+","+last+"\n");
                first=-1; last=-1;
            }
        }
        if (first!=-1) bw.write("chr"+chr+","+first+","+last+"\n");
    }
    
    
    /**
     * Prints statistics on the intersection between regions with >=minValue 
     * overlapping calls and TRF regions.
     */
    private static final void highOverlapIntersectTrf(int chr, int minValue) {
        int i;
        int chrLength;
        double numerator, denominator1, denominator2;
        boolean[] mask;
        int[][] histogram;
        final int SV_TYPE = 4;
        
        if (chr==21) { histogram=histogram_chr21; chrLength=CHR21_LENGTH; mask=trf_mask_chr21; }
        else { histogram=histogram_chr22; chrLength=CHR22_LENGTH; mask=trf_mask_chr22; }
        numerator=0; denominator1=0; denominator2=0;
        for (i=0; i<chrLength; i++) {
            if (mask[i]) denominator1++;
            if (histogram[SV_TYPE][i]>=minValue) denominator2++;
            if (mask[i] && histogram[SV_TYPE][i]>=minValue) numerator++;
        }
        System.err.println("chr"+chr+": nTrfPositions="+denominator1+" nHighOverlapPositions="+denominator2);
        System.err.println("intersection="+numerator+" ("+(100*numerator/denominator1)+"% of nTrfPositions)  ("+(100*numerator/denominator2)+"% of nHighOverlapPositions)");        
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