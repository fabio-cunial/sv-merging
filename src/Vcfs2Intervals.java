import java.io.*;
import java.util.*;
import java.awt.*;
import java.awt.image.*;
import java.awt.geom.*;
import javax.imageio.*;
import java.text.*;
    
    
public class Vcf2Intervals {
    /**
     * INS are represented as intervals of constant length centered at the
     * insertion position.
     */
    private static final int INS_SLACK = 100;  // Arbitrary
    
    private static BufferedImage image;
    private static int N_ROWS, N_COLUMNS;
    private static final int COLOR_BACKGROUND = 0x00F3F3F3;
    private static final int COLOR_BACKGROUND_LIGHT = 0x00FFFFFF;
    private static final int COLOR_BACKGROUND_LINE = 0x00999999;
    private static final int COLOR_REPEAT_SATELLITE = COLOR_BACKGROUND_LINE;
    private static final int COLOR_REPEAT = 0x00CCCCCC;
    private static final int COLOR_CODING = 0x00FF0000;
    private static final int COLOR_SV_START = COLOR_BACKGROUND_LIGHT;
    
    private static final int[] TYPE_COLORS = new int[] {0x00000000, 0x001E86F4, 0x001ABAFA, 0x000C9D46, 0x003BD8BD, 0x00E66703, 0x00F66B0B, 0x008876FE, 0x001D1D1F};
    private static final int QUANTUM = 10;
    private static final int HEIGHT_PIXELS_FILE = 3;
    private static final int HEIGHT_PIXELS_REPEAT = 12;
    
    /**
     * For every track (rows), the number of SVs that cover each position
     * (columns).
     */
    private static int[][] histogram;
    

    public static void main(String[] args) throws Exception {
        final String CHR = args[0];
        final int FROM_POS = Integer.parseInt(args[1]);
        final int TO_POS = Integer.parseInt(args[2]);
        final String REPEAT_MASKER_FILE = args[3];
        final int REPEAT_MASKER_NROWS = Integer.parseInt(args[4]);
        final String TRF_FILE = args[5];
        final int TRF_FILE_NROWS = Integer.parseInt(args[6]);
        final String CODING_FILE = args[7];
        final int CODING_FILE_NROWS = Integer.parseInt(args[8]);
        final String INPUT_FILE_BCFTOOLS_MERGE = args[9];
        final String INPUT_FILE_TRUVARI = args[10];
        final String INPUT_FILE_JASMINE = args[11];
        final String INPUT_FILE_SVIMMER = args[12];
        final String OUTPUT_FILE = args[13];
        
        final boolean ONLY_PASS = false;
        N_COLUMNS=(TO_POS-FROM_POS)/QUANTUM+1;
        
        boolean sameChromosome;
        int i, j, x, y, y0, yMax, yPrime;
        int type, length, color, position;
        String str;
        BufferedReader br;
        int[] maxRows;
        String[] tokens;
        
        
        System.err.println("Estimating number of rows...");
        histogram = new int[4][N_COLUMNS];
        drawVCF(0,INPUT_FILE_SVIMMER,CHR,ONLY_PASS,FROM_POS,TO_POS,0,0);
        drawVCF(1,INPUT_FILE_JASMINE,CHR,ONLY_PASS,FROM_POS,TO_POS,0,0);
        drawVCF(2,INPUT_FILE_TRUVARI,CHR,ONLY_PASS,FROM_POS,TO_POS,0,0);
        drawVCF(3,INPUT_FILE_BCFTOOLS_MERGE,CHR,ONLY_PASS,FROM_POS,TO_POS,0,0);
        maxRows = new int[histogram.length];
        for (i=0; i<histogram.length; i++) {
            maxRows[i]=0;
            for (j=0; j<histogram[i].length; j++) {
                if (histogram[i][j]>maxRows[i]) maxRows[i]=histogram[i][j];
            }
        }
        for (i=0; i<maxRows.length; i++) {
            maxRows[i]*=HEIGHT_PIXELS_FILE+2;
            System.err.println("maxRows["+i+"]="+maxRows[i]);
        }
        histogram=null;
        N_ROWS=0;
        for (i=0; i<maxRows.length; i++) N_ROWS+=maxRows[i]+1;
        N_ROWS+=3*((HEIGHT_PIXELS_REPEAT+2)*10+1);  // 10 arbitrary (for repeats).
        
        System.err.println("Allocating image...");
        image = new BufferedImage(N_COLUMNS,N_ROWS,BufferedImage.TYPE_INT_RGB);
		for (x=0; x<N_COLUMNS; x++) {
		    for (y=0; y<N_ROWS; y++) image.setRGB(x,y,COLOR_BACKGROUND);
		}
        y=0;
        
        System.err.println("Printing: svimmer");
        drawVCF(-1,INPUT_FILE_SVIMMER,CHR,ONLY_PASS,FROM_POS,TO_POS,y,y+maxRows[0]-1);
        y+=maxRows[0];
        for (x=0; x<N_COLUMNS; x++) image.setRGB(x,y,COLOR_BACKGROUND_LINE);
        y++;
        
        System.err.println("Printing: jasmine");
        drawVCF(-1,INPUT_FILE_JASMINE,CHR,ONLY_PASS,FROM_POS,TO_POS,y,y+maxRows[1]-1);
        y+=maxRows[1];
        for (x=0; x<N_COLUMNS; x++) image.setRGB(x,y,COLOR_BACKGROUND_LINE);
        y++;
        
        System.err.println("Printing: truvari");
        drawVCF(-1,INPUT_FILE_TRUVARI,CHR,ONLY_PASS,FROM_POS,TO_POS,y,y+maxRows[2]-1);
        y+=maxRows[2];
        for (x=0; x<N_COLUMNS; x++) image.setRGB(x,y,COLOR_BACKGROUND_LINE);
        y++;
        
        System.err.println("Printing: RepeatMasker annotations");
        yMax=drawRepeatMaskerAnnotations(REPEAT_MASKER_FILE,REPEAT_MASKER_NROWS,CHR,FROM_POS,TO_POS,y,N_ROWS-1);
        y=yMax+1;
        for (x=0; x<N_COLUMNS; x++) image.setRGB(x,y,COLOR_BACKGROUND_LINE);
        y++;
        
        System.err.println("Printing: TRF annotations");
        yMax=drawTrfAnnotations(TRF_FILE,TRF_FILE_NROWS,CHR,FROM_POS,TO_POS,y,N_ROWS-1);
        y=yMax+1;
        for (x=0; x<N_COLUMNS; x++) image.setRGB(x,y,COLOR_BACKGROUND_LINE);
        y++;
        
        System.err.println("Printing: coding regions");
        yMax=drawCodingRegions(-1,CODING_FILE,CHR,FROM_POS,TO_POS,y,N_ROWS-1);
        y=yMax+1;
        for (x=0; x<N_COLUMNS; x++) image.setRGB(x,y,COLOR_BACKGROUND_LINE);
        y++;
        
        System.err.println("Printing: bcftools merge");
        yMax=y+maxRows[3]-1;
        if (yMax>N_ROWS-1) yMax=N_ROWS-1;
        drawVCF(-1,INPUT_FILE_BCFTOOLS_MERGE,CHR,ONLY_PASS,FROM_POS,TO_POS,y,yMax);
        
        ImageIO.write(image,"png",new File(OUTPUT_FILE));
    }
    
    
    private static final int vcf2Intervals(int fileID, String path, String chr, boolean onlyPass, int frameFromX, int frameToX, int frameFromY, int frameToY) throws IOException {
        final int DISTANCE_THRESHOLD = 50;  // Arbitrary
        boolean sameChromosome;
        int x, y, y0, yMax, yPrime;
        int type, length, color, position;
        String str;
        BufferedReader br;
        String[] tokens;
        
        br = new BufferedReader(new FileReader(path));
        str=br.readLine();
        chr=""; first=-1; last=-1;
        while (str!=null) {
            if (str.charAt(0)==COMMENT) {
                str=br.readLine();
                continue;
            }
            tokens=str.split("\t");
            position=Integer.parseInt(tokens[1]);
            if (chr.length()>0 && (!tokens[0].equalsIgnoreCase(chr) || position>last+DISTANCE_THRESHOLD)) {
                
                
                
                
                str=br.readLine();
                continue;
            }
            if (sameChromosome && position>frameToX) break;
			if (onlyPass && !tokens[6].equalsIgnoreCase(PASS_STR)) {
				str=br.readLine();
				continue;
			}
            type=getType_infoField(getField(tokens[7],SVTYPE_STR));     
            if (type!=TYPE_BREAKEND) {
                color=TYPE_COLORS[type-1];
                if (type==TYPE_INSERTION) {
                    length=INS_SLACK*2;
                    position-=INS_SLACK;
                    if (position<0) position=0;
                }
                else {
                    length=Integer.parseInt(getField(tokens[7],SVLEN_STR));
                    if (length<0) length=-length;
                }
                yPrime=drawSV(fileID,frameFromX,frameFromY,frameToY,position,position+length-1,color,HEIGHT_PIXELS_FILE);
                if (yPrime>yMax) yMax=yPrime;
            }
            str=br.readLine();
        }
        br.close();
        return yMax;
    }
    
    
    /**
     * @param path the file is assumed to have no header and rows in the 
     * following format: CHROM, POS, END, ID, SVTYPE, SVLEN. All the SVs in the
     * file are assumed to have passed the filters (FILTER=PASS in a VCF).
     */
    private static final int drawBED(int fileID, String path, String chr, int frameFromX, int frameToX, int frameFromY, int frameToY) throws IOException {
        boolean sameChromosome;
        int x, y, y0, yMax, yPrime;
        int type, length, color, position;
        String str;
        BufferedReader br;
        String[] tokens;
        
        br = new BufferedReader(new FileReader(path));
        str=br.readLine(); yMax=-1;
        while (str!=null) {
            tokens=str.split("\t");
            position=Integer.parseInt(tokens[1]);
            sameChromosome=tokens[0].equalsIgnoreCase(chr);
            if (!sameChromosome || position<frameFromX) {
                str=br.readLine();
                continue;
            }
            if (sameChromosome && position>frameToX) break;
            type=getType_infoField(tokens[4]);     
            if (type!=TYPE_INSERTION && type!=TYPE_BREAKEND) {
                color=TYPE_COLORS[type-1];
                length=Integer.parseInt(tokens[5]);
                if (length<0) length=-length;
                yPrime=drawSV(fileID,frameFromX,frameFromY,frameToY,position,position+length-1,color,HEIGHT_PIXELS_FILE);
                if (yPrime>yMax) yMax=yPrime;
            }
            str=br.readLine();
        }
        br.close();
        if (yMax==-1) yMax=frameFromY;
        return yMax;
    }
    

    /**
     * @param fileID if >=0, the procedure just increments global variable
     * $histogram$;
     * @return the maximum Y that has been drawn, or -1 if $fileID>=0$.
     */
    private static final int drawSV(int fileID, int frameFromX, int frameFromY, int frameToY, int from, int to, int color, int heightPixels) {
        final int MASK = 0x00FFFFFF;
        boolean found;
        int x, y, i, j;
        int fromX, toX;
        
        fromX=(from-frameFromX)/QUANTUM; toX=(to-frameFromX)/QUANTUM;
        if (toX>=N_COLUMNS) toX=N_COLUMNS-1;
        if (histogram!=null) {
            for (x=fromX; x<=toX; x++) histogram[fileID][x]++;
        }
        if (fileID>=0) return -1;
        
        y=frameFromY;
		while (true) {
            if (y>frameToY) return frameToY;
			found=false;
			for (x=fromX; x<=toX; x++) {
                for (j=0; j<heightPixels+2; j++) {
                    if (y+j>frameToY) return frameToY;
    				if ((image.getRGB(x,y+j)&MASK)!=(COLOR_BACKGROUND&MASK)) {
    					found=true;
    					break;
    				}
                }
			}
			if (found) y++;
			else break;
		}
        for (x=fromX; x<=toX; x++) image.setRGB(x,y,COLOR_BACKGROUND_LIGHT);
        for (i=0; i<heightPixels; i++) {
            if (y+1+i>frameToY) return frameToY;
            for (x=fromX; x<=toX; x++) image.setRGB(x,y+1+i,color);
        }
        image.setRGB(fromX,y+1+heightPixels/2,COLOR_SV_START);
        if (y+heightPixels+1<=frameToY) {
            for (x=fromX; x<=toX; x++) image.setRGB(x,y+heightPixels+1,COLOR_BACKGROUND_LIGHT);
            return y+heightPixels+1;
        }
        else return frameToY;
    }

    
	/**
	 * @param path assumed to be a cleaned version of the original RepeatMasker 
	 * file from the UCSC Genome Browser, with runs of spaces replaced by a 
	 * single comma, without header, and already sorted by contig (according to
	 * $VCF2gfa.string2contig()$: this is not the lex order) and starting 
	 * position;
     * @return the largest value of Y used for drawing.
	 */
	private static final int drawRepeatMaskerAnnotations(String path, int nRecords, String chr, int frameFromX, int frameToX, int frameFromY, int frameToY) throws IOException {
        boolean isSatellite, sameChromosome;
		int i, j, k, y;
		int referenceStart, referenceEnd, yMax;
		String str;
		BufferedReader br;
		String[] tokens;

		br = new BufferedReader(new FileReader(path));
		str=br.readLine(); yMax=-1;
		while (str!=null) {
			tokens=str.split(",");
            referenceStart=Integer.parseInt(tokens[5]);
            sameChromosome=tokens[4].equalsIgnoreCase(chr);
            if (!sameChromosome || referenceStart<frameFromX) {
                str=br.readLine();
                continue;
            }
            if (sameChromosome && referenceStart>frameToX) break;
            str=tokens[10].toLowerCase();
			isSatellite=str.indexOf("satellite")>=0 || str.indexOf("sat")>=0 || str.indexOf("simple")>=0 || str.indexOf("low_complexity")>=0;
			referenceEnd=Integer.parseInt(tokens[6]);
            y=drawSV(-1,frameFromX,frameFromY,frameToY,referenceStart,referenceEnd,isSatellite?COLOR_REPEAT_SATELLITE:COLOR_REPEAT,HEIGHT_PIXELS_REPEAT);
            if (y>yMax) yMax=y;
			str=br.readLine();
		}
		br.close();
        if (yMax==-1) yMax=frameFromY;
        return yMax;
    }


    /**
     * @return the largest value of Y used for drawing.
     */
	private static final int drawTrfAnnotations(String path, int nRecords, String chr, int frameFromX, int frameToX, int frameFromY, int frameToY) throws IOException {
        boolean isSatellite, sameChromosome;
		int i, j, k, y;
		int referenceStart, referenceEnd, yMax;
		String str;
		BufferedReader br;
		String[] tokens;

		br = new BufferedReader(new FileReader(path));
		str=br.readLine(); yMax=-1;
		while (str!=null) {
			tokens=str.split(",");
            sameChromosome=tokens[0].equalsIgnoreCase(chr);
            referenceStart=Integer.parseInt(tokens[1]);
            if (!sameChromosome || referenceStart<frameFromX) {
                str=br.readLine();
                continue;
            }
            if (sameChromosome && referenceStart>frameToX) break;
            referenceEnd=Integer.parseInt(tokens[2]);
            y=drawSV(-1,frameFromX,frameFromY,frameToY,referenceStart,referenceEnd,COLOR_REPEAT_SATELLITE,HEIGHT_PIXELS_REPEAT);
            if (y>yMax) yMax=y;
			str=br.readLine();
		}
		br.close();
        if (yMax==-1) yMax=frameFromY;
        return yMax;
    }
    
    
    /**
     * @param path the file is assumed to have one header line and rows in the 
     * following format: CHROM, POS, END.
     */
    private static final int drawCodingRegions(int fileID, String path, String chr, int frameFromX, int frameToX, int frameFromY, int frameToY) throws IOException {
        boolean sameChromosome;
        int x, y, y0, yMax, yPrime;
        int type, length, color, position;
        String str;
        BufferedReader br;
        String[] tokens;
        
        br = new BufferedReader(new FileReader(path));
        str=br.readLine();  // Skipping header
        str=br.readLine(); yMax=-1;
        while (str!=null) {
            tokens=str.split("\t");
            position=Integer.parseInt(tokens[1]);
            sameChromosome=tokens[0].equalsIgnoreCase(chr);
            if (!sameChromosome || position<frameFromX) {
                str=br.readLine();
                continue;
            }
            if (sameChromosome && position>frameToX) break;
            yPrime=drawSV(fileID,frameFromX,frameFromY,frameToY,position,Integer.parseInt(tokens[2]),COLOR_CODING,HEIGHT_PIXELS_REPEAT);
            if (yPrime>yMax) yMax=yPrime;
            str=br.readLine();
        }
        br.close();
        if (yMax==-1) yMax=frameFromY;
        return yMax;
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
    
    
	/**
	 * @return NULL if $field$ does not occur in $str$.
	 */
	public static final String getField(String str, String field) {
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
	 * @return -1 iff the type cannot be determined.
	 */
	private static final int getType_infoField(String type) {
		if (type==null || type.length()==0) return -1;
		if ( type.equalsIgnoreCase(DEL_STR) || 
			 type.equalsIgnoreCase(DEL_ME_STR)
		   ) return TYPE_DELETION;
		else if (type.equalsIgnoreCase(DEL_INV_STR)) return TYPE_DEL_INV;
		else if ( type.equalsIgnoreCase(INS_STR) || 
			      type.equalsIgnoreCase(INS_ME_STR) || 
				  type.equalsIgnoreCase(INS_NOVEL_STR)
				) return TYPE_INSERTION;
		else if ( type.equalsIgnoreCase(DUP_STR) ||
			      type.equalsIgnoreCase(DUP_TANDEM_STR) ||
				  type.equalsIgnoreCase(DUP_INT_STR)
			    ) return TYPE_DUPLICATION;
		else if (type.equalsIgnoreCase(INV_STR)) return TYPE_INVERSION;
		else if (type.equalsIgnoreCase(INV_DUP_STR)) return TYPE_INV_DUP;
		else if (type.equalsIgnoreCase(CNV_STR)) return TYPE_CNV;
		else if (type.equalsIgnoreCase(BND_STR)) return TYPE_BREAKEND;
		else if (type.equalsIgnoreCase(TRA_STR)) return TYPE_TRANSLOCATION;
		else return -1;
	}

}