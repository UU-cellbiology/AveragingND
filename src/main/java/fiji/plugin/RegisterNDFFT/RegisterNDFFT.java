package fiji.plugin.RegisterNDFFT;


import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import io.scif.img.ImgIOException;
import io.scif.img.ImgOpener;
import net.imglib2.exception.IncompatibleTypeException;
import net.imglib2.img.ImagePlusAdapter;
import net.imglib2.img.Img;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;


/**
 * Perform template matching by finding maximum of normalized cross-correlation
 * using convolution in the Fourier domain.
 *
 * @author Eugene Katrukha
 * 
 */
public class RegisterNDFFT implements PlugIn 
{
	
	public static int defaultImg1 = 0;
	public static int defaultImg2 = 1;
	public static double dMaxFraction = 0.5;
	public static boolean bShowCC = true;
	public static boolean bRegisterTemplate = true;
	public long [] finShift;
	public int regChannel1 =0;
	public int regChannel2 =0;
	public boolean multiCh = false;

	@Override
	public void run(String arg) {
		
		int i;

		//this part is honestly stolen from "Pairwise Stitching" plugin
		//https://github.com/fiji/Stitching/blob/master/src/main/java/plugin/Stitching_Pairwise.java
		// get list of image stacks
		final int[] idList = WindowManager.getIDList();		

		if ( idList == null || idList.length < 2 )
		{
			IJ.error( "You need at least two open images." );
			return;
		}

		final String[] imgList = new String[ idList.length ];
		for ( i = 0; i < idList.length; ++i )
			imgList[ i ] = WindowManager.getImage(idList[i]).getTitle();
		
		/**
		 * Dialog for choosing the images
		 */
		final GenericDialog gd1 = new GenericDialog( "Images/volumes registration" );
		
		if ( defaultImg1 >= imgList.length || defaultImg2 >= imgList.length )
		{
			defaultImg1 = 0;
			defaultImg2 = 1;
		}
		
		gd1.addChoice("First_image (reference)", imgList, imgList[ defaultImg1 ] );
		gd1.addChoice("Second_image (to register)", imgList, imgList[ defaultImg2 ] );
		gd1.addNumericField("Maximum shift (fraction, 0-1 range)", Prefs.get("RegisterNDFFT.dMaxFraction", 0.5), 3);
		gd1.addCheckbox("Show cross-correlation", Prefs.get("RegisterNDFFT.bShowCC", false));
		gd1.addCheckbox("Register template", Prefs.get("RegisterNDFFT.bRegisterTemplate", false));
				
		gd1.showDialog();
		
		if ( gd1.wasCanceled() )
			return;
		
		ImagePlus imp1 = WindowManager.getImage( idList[ defaultImg1 = gd1.getNextChoiceIndex() ] );		
		ImagePlus imp2 = WindowManager.getImage( idList[ defaultImg2 = gd1.getNextChoiceIndex() ] );	
		
		
		dMaxFraction  = gd1.getNextNumber();
		Prefs.set("RegisterNDFFT.dMaxFraction", dMaxFraction);
		bShowCC  = gd1.getNextBoolean();
		Prefs.set("RegisterNDFFT.bShowCC", bShowCC);
		bRegisterTemplate  = gd1.getNextBoolean();
		Prefs.set("RegisterNDFFT.bRegisterTemplate", bRegisterTemplate);
		
		// create channel selector
		final int numChannels1 = imp1.getNChannels();
		final int numChannels2 = imp2.getNChannels();
		
		final String[] channels1 = new String[ numChannels1 ];
		final String[] channels2 = new String[ numChannels2 ];
		if(numChannels1>1 && numChannels2>1)
		{
			multiCh = true;
		}
		
		
		if(multiCh)
		{
			for ( int c = 0; c < channels1.length; ++c )
				channels1[ c ] = "Use channel " + Integer.toString(c+1);
			for ( int c = 0; c < channels2.length; ++c )
				channels2[ c ] = "Use channel " + Integer.toString(c+1);
			
			final GenericDialog gd2 = new GenericDialog( "Register multichannel" );

			gd2.addChoice( "Reference_image_channel", channels1, channels1[ 0 ] );
			gd2.addChoice( "Registered_image_channel", channels2, channels2[ 0 ] );
			gd2.showDialog();

			if ( gd2.wasCanceled() )
				return;
			regChannel1 = gd2.getNextChoiceIndex();
			regChannel2 = gd2.getNextChoiceIndex();
			
		}
		
		final Img< FloatType > image_in = ImagePlusAdapter.convertFloat(imp1);
		final Img< FloatType > template_in = ImagePlusAdapter.convertFloat(imp2);
		GenNormCC normCC = new GenNormCC();
		boolean bNormCCcalc=false;
		if(multiCh)
		{	
			bNormCCcalc= normCC.caclulateGenNormCC(Views.hyperSlice(image_in, 2, regChannel1), Views.hyperSlice(template_in, 2, regChannel2), dMaxFraction , bShowCC);//, bRegisterTemplate);
			//finData=GenNormCC.caclulateGenNormCC(Views.hyperSlice(image_in, 2, regChannel1), Views.hyperSlice(template_in, 2, regChannel2), dMaxFraction , bShowCC);//, bRegisterTemplate);
		}
		else					
		{
			bNormCCcalc = normCC.caclulateGenNormCC(image_in, template_in, dMaxFraction , bShowCC);//, bRegisterTemplate);	
			//finData=GenNormCC.caclulateGenNormCC(image_in, template_in, dMaxFraction , bShowCC);//, bRegisterTemplate);
		}
		if(!bNormCCcalc)
		{
			return;
		}
		
		finShift= normCC.dShift;
		

		ResultsTable ptable = ResultsTable.getResultsTable();
		//RegisterTranslation reg = new RegisterTranslation(image_in.numDimensions());
		
		//reg.registerTranslation(image_in, template_in);
		//ptable_lock.lock();
		ptable.incrementCounter();
		ptable.addValue("norm_CC_coeff", normCC.dMaxCC);
		for(i=0;i<finShift.length;i++)
		{
			ptable.addValue("shift_coord_"+Integer.toString(i),finShift[i]);	
		}
		ptable.show("Results");
		
		if(bRegisterTemplate)
		{
				registerShowTemplate(finShift,template_in,image_in, imp2.getTitle(),imp2.getCalibration(),multiCh);
		}
		//ptable_lock.unlock();
	}
	
	private void registerShowTemplate(final long[] finData, final Img<FloatType> template, final Img<FloatType> image, final String sTitle, final Calibration cal, final boolean bCh) 
	{
		int nDim = template.numDimensions();
		long [] shift = new long [nDim];
		int i,j;
		
		//"jumping" over color channel, since it is xyczt and our shift is xyz
		j=0;
		for (i=0;i<nDim;i++)
		{
			shift[i]=finData[j];
			if(bCh && i==1)
			{
				shift[2]=0;
				i++;
			}
			j++;
		}
		

		ImagePlus registeredIP = MiscUtils.wrapFloatImgCal(Views.interval(Views.translate(Views.extendZero(template), shift),image), "registered_"+sTitle, cal, bCh);
		
		registeredIP.show();
	}

	public static void main( final String[] args ) throws ImgIOException, IncompatibleTypeException
	{
		// open an ImageJ window
		new ImageJ();
		
		// open with SCIFIO ImgOpener as FloatTypes
		ImgOpener io = new ImgOpener();
		
		final Img< FloatType > image_in = io.openImgs( "linetest2Dn_ve.tif",
			new FloatType() ).get( 0 );
		//final Img< FloatType > template_in = io.openImgs( "Zstack 1.1-1-1.tif",
		final Img< FloatType > template_in = io.openImgs( "linetest2Dn_crop.tif",
			new FloatType() ).get( 0 );
		
			
		
		/*
		final Img< FloatType > image_in = io.openImgs( "s001.tif",
				new FloatType() ).get( 0 );
		final Img< FloatType > template_in = io.openImgs( "s002.tif",
				new FloatType() ).get( 0 );
		*/
		long[] imdim= new long[image_in.numDimensions()];
		image_in.dimensions(imdim);
		

/*	final Img< FloatType > image_in = io.openImgs( "linetest2Dn.tif",
				new FloatType() ).get( 0 );
		final Img< FloatType > template_in = io.openImgs( "linetest2Dn_shift.tif",
				new FloatType() ).get( 0 );
*/
		GenNormCC normCC = new GenNormCC();
		normCC.caclulateGenNormCC(image_in, template_in, 0.50, true);//, true);
		//GenNormCC.caclulateGenNormCC(Views.hyperSlice(image_in, 2, 2), Views.hyperSlice(template_in, 2, 2), 0.4, true, true);
		//RegisterTranslation reg = new RegisterTranslation(image_in .numDimensions());
		// run the example
		//reg.registerTranslation(image_in, template_in);

		
	
	}
	//void showRegisteredTemplate()

}
