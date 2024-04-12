package registerNDFFT;


import java.awt.AWTEvent;
import java.awt.Checkbox;
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
public class RegisterNDFFT implements PlugIn, DialogListener
{
	
	public static int defaultImg1 = 0;
	public static int defaultImg2 = 1;
	public static boolean bShowCC = true;
	public static boolean bRegisterTemplate = true;
	public long [] finShift;
	public int regChannel1 = 0;
	public int regChannel2 = 0;
	public boolean multiCh = false;
	public boolean bExcludeZeros = false;
	public int nConstrainReg = 0;
	Label [] limName;
	TextField [] limVal;
	int nDimReg;

	@Override
	public void run(String arg) {
		
		int i;
		
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
		for ( i = 0; i < idList.length; ++i )
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
		
		String sDims = MiscUtils.getDimensionsText(imp1);
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
		gd1.addCheckbox("Exclude zero values?", Prefs.get("RegisterNDFFT.bExcludeZeros", false));		
		gd1.addCheckbox("Show cross-correlation?", Prefs.get("RegisterNDFFT.bShowCC", false));
		gd1.addCheckbox("Register template?", Prefs.get("RegisterNDFFT.bRegisterTemplate", false));
		String sCurrChoice = Prefs.get("RegisterNDFFT.sConstrain", "No");
		gd1.addChoice("Constrain registration?", limitsReg, sCurrChoice);
		
		for (i=0;i<nDimReg;i++)
		{
			switch (sCurrChoice)
			{
				case "No":
					gd1.addNumericField("No max "+sDims.charAt(i)+" limit", 0.0, 3);
					break;
				case "by voxels":
					gd1.addNumericField(sDims.charAt(i)+" limit (px)", Prefs.get("RegisterNDFFT.dMax"+sDims.charAt(i)+"px", 10.0), 3);
					break;
				case "by image fraction":
					gd1.addNumericField(sDims.charAt(i)+" limit (0-1)", Prefs.get("RegisterNDFFT.dMax"+sDims.charAt(i)+"fr", 0.5), 3);
					break;
					
			}
			limName[i] = gd1.getLabel();
			limVal[i] = (TextField)gd1.getNumericFields().get(0);	
			if(sCurrChoice.equals("No"))
			{
				limVal[i].setEnabled(false);
			}
		}


		gd1.addDialogListener(this);
		gd1.showDialog();
		
		if ( gd1.wasCanceled() )
			return;
		
		bExcludeZeros  = gd1.getNextBoolean();
		Prefs.set("RegisterNDFFT.bExcludeZeros", bExcludeZeros);
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
				for(i=0;i<nDimReg;i++)
				{
					dLimits[i]=Math.abs(gd1.getNextNumber());
				}
				Prefs.set("RegisterNDFFT.dMax"+sDims.charAt(i)+"fr",dLimits[i]);
			}
			else
			{
				for(i=0;i<nDimReg;i++)
				{
					dLimits[i]=Math.min(Math.abs(gd1.getNextNumber()), 1.0);
				}
				Prefs.set("RegisterNDFFT.dMax"+sDims.charAt(i)+"px",dLimits[i]);
			

			}
		}
		double [] lim_fractions = null;
		FinalInterval limInterval = null;
		if(nConstrainReg == 1)
		{
			long[] minI;
			long[] maxI;
			if(bZPresent)
			{
				minI = new long[3];
				maxI = new long[3];
				minI[2] = (long) ((-1.0)*dLimZ);
				maxI[2] = (long) (dLimZ);
			}
			else
			{
				minI = new long[2];
				maxI = new long[2];
			}
			maxI[0] = (long) dLimX;
			maxI[1] = (long) dLimY;
			minI[0] = (long) ((-1.0)*dLimX);
			minI[1] = (long) ((-1.0)*dLimY);
			limInterval = new FinalInterval(minI, maxI);
		}
		if(nConstrainReg == 2)
		{
			if(bZPresent)
			{
				lim_fractions = new double[3];
				lim_fractions[2] = dLimZ;
			}
			else
			{
				lim_fractions = new double[2];
			}
			lim_fractions[0] = dLimX;
			lim_fractions[1] = dLimY;
			
		}

		
		final Img< FloatType > image_in = ImagePlusAdapter.convertFloat(imp1);
		final Img< FloatType > template_in = ImagePlusAdapter.convertFloat(imp2);
		GenNormCC normCC = new GenNormCC();
		normCC.bExcludeZeros = bExcludeZeros;
		
		normCC.lim_fractions = lim_fractions;
		normCC.limInterval = limInterval;
		
		boolean bNormCCcalc=false;
		
		if(multiCh)
		{			
			bNormCCcalc= normCC.caclulateGenNormCC(Views.hyperSlice(image_in, 2, regChannel1), Views.hyperSlice(template_in, 2, regChannel2), bShowCC);//, bRegisterTemplate);
			
			//finData=GenNormCC.caclulateGenNormCC(Views.hyperSlice(image_in, 2, regChannel1), Views.hyperSlice(template_in, 2, regChannel2), dMaxFraction , bShowCC);//, bRegisterTemplate);
		}
		else					
		{
			bNormCCcalc = normCC.caclulateGenNormCC(image_in, template_in, bShowCC);//, bRegisterTemplate);	
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
	@Override
	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
		
		if(e!=null)
		{
			DecimalFormatSymbols symbols = new DecimalFormatSymbols();
			symbols.setDecimalSeparator('.');
			DecimalFormat df1 = new DecimalFormat ("#.##", symbols);
			Choice limit = (Choice) gd.getChoices().get(0);
			if(e.getSource()==limit)
			{
				switch (limit.getSelectedIndex())
				{
					case 0:
						xName.setText("No X limit");										
						yName.setText("No Y limit");
						xVal.setEnabled(false);
						yVal.setEnabled(false);
						if(bZPresent)
						{
							zName.setText("No Z limit");
							zVal.setEnabled(false);
						}
						break;
					case 1:
						xName.setText("X limit (px)");
						yName.setText("Y limit (px)");
						xVal.setEnabled(true);
						yVal.setEnabled(true);
						xVal.setText(df1.format(Prefs.get("RegisterNDFFT.dMaxXpx", 10.0)));
						yVal.setText(df1.format(Prefs.get("RegisterNDFFT.dMaxYpx", 10.0)));
						if(bZPresent)
						{
							zName.setText("Z limit (px)");
							
							zVal.setEnabled(true);
							zVal.setText(df1.format(Prefs.get("RegisterNDFFT.dMaxZpx", 10.0)));
						}
						break;
					case 2:
						xName.setText("X limit (0-1)");
						yName.setText("Y limit (0-1)");
						xVal.setEnabled(true);
						yVal.setEnabled(true);
						xVal.setText(df1.format( Prefs.get("RegisterNDFFT.dMaxXfr", 0.5)));
						yVal.setText(df1.format( Prefs.get("RegisterNDFFT.dMaxYfr", 0.5)));
						if(bZPresent)
						{
							zName.setText("Z limit (0-1)");
							zVal.setEnabled(true);
							zVal.setText(df1.format( Prefs.get("RegisterNDFFT.dMaxZfr", 0.5)));
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
		
		IJ.open("/home/eugene/Desktop/projects/RegisterNDFFT/4d/HyperStack.tif");
		IJ.open("/home/eugene/Desktop/projects/RegisterNDFFT/4d/HyperStack-1.tif");
		
		RegisterNDFFT test = new RegisterNDFFT();
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
