package averagingND;


import java.awt.AWTEvent;
import java.awt.Label;
import java.awt.TextField;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.awt.Choice;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import io.scif.img.ImgIOException;

import net.imglib2.FinalInterval;
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
public class RegisterSingleND implements PlugIn, DialogListener
{
	
	public static int defaultImg1 = 0;
	public static int defaultImg2 = 1;
	public static boolean bShowCC = true;
	public static boolean bRegisterTemplate = true;
	public long [] finShift;
	public int regChannel1 = 0;
	public int regChannel2 = 0;
	public boolean multiCh = false;
	public boolean bZeroMask = false;
	public int nConstrainReg = 0;
	
	/** choice UI for constrain type **/
	Choice limitCh;
	
	/** labels of constrain axes **/
	Label [] limName;
	
	/** values of constrain axes **/
	TextField [] limVal;
	
	int nDimReg;
	String sDims;
	boolean bCenteredLimit = false;

	@Override
	public void run(String arg) {
		
		int d;
		
		double [] dLimits;

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
		for ( int i = 0; i < idList.length; ++i )
			imgList[ i ] = WindowManager.getImage(idList[i]).getTitle();
		
		/**
		 * Dialog for choosing the images
		 */
		final GenericDialog gdImages = new GenericDialog( "Images/volumes registration" );
	
		if ( defaultImg1 >= imgList.length || defaultImg2 >= imgList.length )
		{
			defaultImg1 = 0;
			defaultImg2 = 1;
		}
		
		gdImages.addChoice("First_image (reference)", imgList, imgList[ defaultImg1 ] );
		gdImages.addChoice("Second_image (template)", imgList, imgList[ defaultImg2 ] );
		
		gdImages.showDialog();
				
		if ( gdImages.wasCanceled() )
			return;
		
		//TODO here I should check that they have the same dimensions
		
		ImagePlus imp1 = WindowManager.getImage( idList[ defaultImg1 = gdImages.getNextChoiceIndex() ] );		
		ImagePlus imp2 = WindowManager.getImage( idList[ defaultImg2 = gdImages.getNextChoiceIndex() ] );	
		
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
			
			final GenericDialog gdCh = new GenericDialog( "Choose registration channel" );

			gdCh.addChoice( "Reference_image_channel", channels1, channels1[ 0 ] );
			gdCh.addChoice( "Template_image_channel", channels2, channels2[ 0 ] );
			gdCh.showDialog();

			if ( gdCh.wasCanceled() )
				return;
			regChannel1 = gdCh.getNextChoiceIndex();
			regChannel2 = gdCh.getNextChoiceIndex();
			
		}
		
		sDims = MiscUtils.getDimensionsText(imp1);
		nDimReg = sDims.length();
		if(multiCh)
		{
			nDimReg--; //remove the C component
		}
		limName = new Label[nDimReg];
		limVal = new TextField[nDimReg];
		dLimits = new double [nDimReg];
		
		final String[] limitsReg = new String[  ] {"No","by voxels", "by image fraction"};
		final GenericDialog gd1 = new GenericDialog( "Registration parameters" );	
		gd1.addCheckbox("Use zero masked CC?", Prefs.get("RegisterNDFFT.bExcludeZeros", false));		
		gd1.addCheckbox("Show cross-correlation?", Prefs.get("RegisterNDFFT.bShowCC", false));
		gd1.addCheckbox("Register template?", Prefs.get("RegisterNDFFT.bRegisterTemplate", false));
		String sCurrChoice = Prefs.get("RegisterNDFFT.sConstrain", "No");
		gd1.addChoice("Constrain registration?", limitsReg, sCurrChoice);
		limitCh = (Choice) gd1.getChoices().lastElement();
		
		for (d=0;d<nDimReg;d++)
		{
			switch (sCurrChoice)
			{
				case "No":
					gd1.addNumericField("No max "+sDims.charAt(d)+" limit", 0.0, 3);
					break;
				case "by voxels":
					gd1.addNumericField(sDims.charAt(d)+" limit (px)", Prefs.get("RegisterNDFFT.dMax"+sDims.charAt(d)+"px", 10.0), 3);
					break;
				case "by image fraction":
					gd1.addNumericField(sDims.charAt(d)+" limit (0-1)", Prefs.get("RegisterNDFFT.dMax"+sDims.charAt(d)+"fr", 0.5), 3);
					break;
					
			}
			limName[d] = gd1.getLabel();
			limVal[d] = (TextField)gd1.getNumericFields().get(d);	
			if(sCurrChoice.equals("No"))
			{
				limVal[d].setEnabled(false);
			}
		}
		gd1.addCheckbox("Image centered constrains?", Prefs.get("RegisterNDFFT.bCenteredLimit", false));

		gd1.addDialogListener(this);
		gd1.showDialog();
		
		if ( gd1.wasCanceled() )
			return;
		
		bZeroMask  = gd1.getNextBoolean();
		Prefs.set("RegisterNDFFT.bExcludeZeros", bZeroMask);
		bShowCC  = gd1.getNextBoolean();
		Prefs.set("RegisterNDFFT.bShowCC", bShowCC);
		bRegisterTemplate  = gd1.getNextBoolean();
		Prefs.set("RegisterNDFFT.bRegisterTemplate", bRegisterTemplate);
		nConstrainReg = gd1.getNextChoiceIndex();
		Prefs.set("RegisterNDFFT.sConstrain", limitsReg[nConstrainReg]);
		if(nConstrainReg!=0)
		{
			if(nConstrainReg == 1)
			{
				for(d=0;d<nDimReg;d++)
				{
					dLimits[d]=Math.abs(gd1.getNextNumber());
					Prefs.set("RegisterNDFFT.dMax"+sDims.charAt(d)+"px",dLimits[d]);
				}
				
			}
			else
			{
				for(d=0;d<nDimReg;d++)
				{
					dLimits[d]=Math.min(Math.abs(gd1.getNextNumber()), 1.0);
					Prefs.set("RegisterNDFFT.dMax"+sDims.charAt(d)+"fr",dLimits[d]);
				}
			}
		}
		bCenteredLimit  = gd1.getNextBoolean();
		Prefs.set("RegisterNDFFT.bCenteredLimit", bCenteredLimit);
		
		double [] lim_fractions = null;
		FinalInterval limInterval = null;
		if(nConstrainReg == 1)
		{
			long[] minI = new long [nDimReg];
			long[] maxI = new long [nDimReg];
			for(d=0;d<nDimReg;d++)
			{
				maxI[d] = (long) dLimits[d];
				minI[d] = (long) ((-1.0)*dLimits[d]);
			}
			limInterval = new FinalInterval(minI, maxI);
		}
		if(nConstrainReg == 2)
		{
			lim_fractions = new double [nDimReg];
			for(d=0;d<nDimReg;d++)
			{
				lim_fractions[d] = dLimits[d];
			}
		}
		
		//convert to RAI
		final Img< FloatType > image_in = ImagePlusAdapter.convertFloat(imp1);
		final Img< FloatType > template_in = ImagePlusAdapter.convertFloat(imp2);
		MaskedNormCC normCC = new MaskedNormCC();
		normCC.bZeroMask = bZeroMask;
		
		normCC.lim_fractions = lim_fractions;
		normCC.limInterval = limInterval;
		normCC.bCenteredLimit = bCenteredLimit;
		
		boolean bNormCCcalc = false;
		
		if(multiCh)
		{			
			bNormCCcalc= normCC.caclulateMaskedNormCC(Views.hyperSlice(image_in, 2, regChannel1), Views.hyperSlice(template_in, 2, regChannel2), bShowCC);//, bRegisterTemplate);
		}
		else					
		{
			bNormCCcalc = normCC.caclulateMaskedNormCC(image_in, template_in, bShowCC);//, bRegisterTemplate);	
		}
		if(!bNormCCcalc)
		{
			return;
		}
		
		finShift = normCC.dShift;
		

		ResultsTable ptable = ResultsTable.getResultsTable();
		ptable.incrementCounter();
		ptable.addValue("norm_CC_coeff", normCC.dMaxCC);
		for(d=0;d<finShift.length;d++)
		{
			ptable.addValue("shift_coord_"+Integer.toString(d),finShift[d]);	
		}
		ptable.show("Results");
		
		if(bRegisterTemplate)
		{
				registerShowTemplate(finShift,template_in,image_in, imp2.getTitle(),imp2.getCalibration(),multiCh);
		}

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
	
	
	@Override
	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
		
		int d;
		
		if(e!=null)
		{
			DecimalFormatSymbols symbols = new DecimalFormatSymbols();
			symbols.setDecimalSeparator('.');
			DecimalFormat df1 = new DecimalFormat ("#.##", symbols);
			
			if(e.getSource()==limitCh)
			{
				switch (limitCh.getSelectedIndex())
				{
					case 0:
						for(d=0;d<nDimReg;d++)
						{
							limName[d].setText("No "+sDims.charAt(d)+" limit");
							limVal[d].setEnabled(false);
						}
						break;
					case 1:
						for(d=0;d<nDimReg;d++)
						{
							limName[d].setText(sDims.charAt(d)+" limit (px)");
							limVal[d].setEnabled(true);
							limVal[d].setText(df1.format(Prefs.get("RegisterNDFFT.dMax"+sDims.charAt(d)+"px", 10.0)));
						}
						break;
					case 2:
						for(d=0;d<nDimReg;d++)
						{
							limName[d].setText(sDims.charAt(d)+" limit (0-1)");
							limVal[d].setEnabled(true);
							limVal[d].setText(df1.format(Prefs.get("RegisterNDFFT.dMax"+sDims.charAt(d)+"fr", 0.5)));

						}
						break;
				}
			}
		}
		return true;
	}

	public static void main( final String[] args ) throws ImgIOException, IncompatibleTypeException
	{
		// open an ImageJ window
		new ImageJ();
		//IJ.open("/home/eugene/Desktop/projects/RegisterNDFFT/single/MAX_089-1.tif");
		//IJ.open("/home/eugene/Desktop/projects/RegisterNDFFT/single/MAX_098-1.tif");
		//IJ.open("/home/eugene/Desktop/projects/RegisterNDFFT/single/MAX_089-1.tif");
		
		//IJ.open("/home/eugene/Desktop/projects/RegisterNDFFT/4d/HyperStack.tif");
		//IJ.open("/home/eugene/Desktop/projects/RegisterNDFFT/4d/HyperStack-1.tif");
		//IJ.open("/home/eugene/Desktop/projects/AveragingND/center/full.tif");
		//IJ.open("/home/eugene/Desktop/projects/AveragingND/center/center.tif");
		IJ.open("/home/eugene/Desktop/projects/AveragingND/centered/full.tif");
		IJ.open("/home/eugene/Desktop/projects/AveragingND/centered/center.tif");

		
		RegisterSingleND test = new RegisterSingleND();
		test.run(null);
		
		// open with SCIFIO ImgOpener as FloatTypes
		//ImgOpener io = new ImgOpener();
		

	/*
		final Img< FloatType > image_in = io.openImgs( "bb1smEC.tif",
			new FloatType() ).get( 0 );

		final Img< FloatType > template_in = io.openImgs( "bb1sm20rotEC.tif",
				new FloatType() ).get( 0 );
			*/	
		
		/*
		final Img< FloatType > image_in = io.openImgs( "wave0_20.tif",
			new FloatType() ).get( 0 );
		final Img< FloatType > template_in = io.openImgs( "wave45_20.tif",
			new FloatType() ).get( 0 );
		*/
		/*
		final Img< FloatType > image_in = io.openImgs( "h1.tif",
				new FloatType() ).get( 0 );
			final Img< FloatType > template_in = io.openImgs( "h2.tif",
				new FloatType() ).get( 0 );
			
			*/
		/*
		final Img< FloatType > image_in = io.openImgs( "polarOriginal.tif",
			new FloatType() ).get( 0 );
		final Img< FloatType > template_in = io.openImgs( "polarOriginal_rot45.tif",
			new FloatType() ).get( 0 );
			*/
		
		 //ImageJFunctions.show(image_in).setTitle("image");
		 //ImageJFunctions.show(template_in).setTitle("template");
		
	/*
		long[] imdim= new long[image_in.numDimensions()];
		image_in.dimensions(imdim);
		

		GenNormCC normCC = new GenNormCC();
		normCC.caclulateGenNormCC(image_in, template_in, 0.50, true);//, true);
*/
		/* 
		RotationCC rotCC =new RotationCC();
		
		rotCC.caclulateRotationFFTCC(image_in, template_in);
		*/
		
	
	}




}
