package averagingND;


import java.awt.AWTEvent;
import java.awt.Choice;
import java.awt.Label;
import java.awt.TextField;
import java.io.File;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Arrays;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.io.DirectoryChooser;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;

import net.imglib2.FinalInterval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;

public class IterativeAveraging implements PlugIn, DialogListener {

	
	/** set of images for averaging and information about them **/
	ImageSet imageSet;
	
	/** shifted images for average calculation**/
	ArrayList<RandomAccessibleInterval< FloatType >> imgs_shift;
	
	/** ND shift of each image **/
	ArrayList<long []> shifts;
	
	/** function holding information about current template **/
	TemplateAveraging template = null;
	
	/** source of images 
	 * 0 - images currently open in ImageJ 
	 * 1 - tif files on disk **/
	public int nInputType = 0;
	
	/** initial template (alignment of images)
	 * 0 - centered images
	 * 1 - zero, origin of coordinates
	 * **/
	public int nIniTemplate = 0;
	
	/** total number of iterations **/
	public int nIterN = 0;

	/** whether or not generate intermediate average **/
	boolean bSaveIntermediate = false;
	
	String sPathGeneralOutput = "";
	
	String sPathIntermediate = "";

	/** whether to use zero masked CC **/
	public boolean bZeroMask = false;
	
	public boolean bOutputInput = false;
	
	String sPathRegistered = "";
	
	/** whether to use simple average or zero masked **/
	public boolean bIgnoreZeroInAveraging = true;

	/** constrains during the averaging
	 * 0 - no constrains
	 * 1 - constrains in pixels
	 * 2 - constrains as a fraction of max displacement 
	 * **/
	public int nConstrainReg = 0;
	
	double [] lim_fractions = null;
	
	FinalInterval limInterval = null;
	
	/** choice UI for constrain type **/
	Choice limitCh;
	
	/** labels of constrain axes **/
	Label [] limName;
	
	/** values of constrain axes **/
	TextField [] limVal;
	
	/** dimensions of dataset for averaging (always 1 channel) **/
	int nDimReg;
	
	/** format of the input dataset XYZTC **/
	String sDims;
	
	/** type of averaging aim
	 *  0 - average
	 *  1 - median
	 *  2 - zero masked average **/
	int nAveragingAim;
	
	boolean bShowSD = false;
	
	boolean bLoadCached = false;
	
	@Override
	public void run(String paramString) {
	
		int i,j,k, iter;
		int d;
		
		ArrayList<RandomAccessibleInterval< FloatType >> imgs_avrg_out;
		
		//double format formatting tool
		DecimalFormatSymbols symbols = new DecimalFormatSymbols();
		symbols.setDecimalSeparator('.');
		DecimalFormat df = new DecimalFormat ("#.########", symbols);
		DecimalFormat df1 = new DecimalFormat ("#.#", symbols);
		
		final String[] sInput = new String[3];
		sInput[0] = "All currently open images";
		sInput[1] = "Tif files (high memory, fast)";
		sInput[2] = "Tif files (low memory, slow)";
		
		final GenericDialog gdFiles = new GenericDialog( "Iterative averaging" );

		gdFiles.addChoice( "Input images:", sInput, Prefs.get("RegisterNDFFT.IA.nInput", sInput[0]) );
		gdFiles.showDialog();
		
		if ( gdFiles.wasCanceled() )
			return;			
		IJ.log("Iterative ND averaging plugin, version " + ConstantsAveragingND.sVersion);

		nInputType = gdFiles.getNextChoiceIndex();
		Prefs.set("RegisterNDFFT.IA.nInput", sInput[nInputType]);
		
		imageSet = new ImageSet();
		
		//init arrays			
		imgs_shift = new ArrayList<RandomAccessibleInterval< FloatType >>();
		shifts = new ArrayList<long []>();
		imageSet.nSource = nInputType;
		
		if(nInputType == ImageSet.OPENED)
		{
			if(!imageSet.initializeFromOpenWindows())
				return;	 
			IJ.log("Using images already opened in ImageJ/Fiji.");
		}
		else
		{
			if(nInputType == ImageSet.LOAD_MEMORY)
			{
				IJ.log("Loading all images to memory.");
			}
			else
			{
				IJ.log("Using cached reading of images.");				
			}

			DirectoryChooser dc = new DirectoryChooser ( "Choose a folder with images.." );
			String sPath = dc.getDirectory();
			if(!imageSet.initializeFromDisk(sPath, ".tif"))
				return;
		}
		
		if(!dialogSettings())
			return;
		
		if(!imageSet.loadAllImages())
			return;
		
		if(nIterN>0)
		{
			IJ.log("Running for "+Integer.toString(nIterN)+" iterations.");
		}		
		final int nImageN = imageSet.nImageN;
		
		//initial shifts depending on template choice
		shifts = initShifts(imageSet.imgs, nIniTemplate);
		
		//create a new shifted array of images
		buildShiftedIntervals(imageSet.imgs, imgs_shift,shifts);
		
		
		template = new TemplateAveraging(nAveragingAim);

		//ArrayList<IntervalView<FloatType>> sumAndCount = null;
		//sumAndCount = AverageWithoutZero.sumAndCountArray(imgs_shift);

		IntervalView<FloatType> currAverageImg;
		if(bSaveIntermediate && nIterN>0)
		{		
			processIntermediate(0);
		}
		
		GenNormCC normCC = new GenNormCC();
		normCC.bVerbose = false;
		normCC.bZeroMask = bZeroMask;
		normCC.lim_fractions = lim_fractions;
		normCC.limInterval = limInterval;
		normCC.bCenteredLimit = true;
		
		ResultsTable ptable = ResultsTable.getResultsTable();
		ptable.reset();
		ResultsTable ptableCC = new ResultsTable();
		ptableCC.setPrecision(8);
		
		double avrgCC;
		double cumShift=0.0;
		
		double oldCumShift = -1.0;
		double oldAvrgCC = Double.NaN;
		
		boolean bConverged = false;
		IJ.showStatus("Iterative averaging...");
		IJ.showProgress(0, (nIterN)*(imageSet.nImageN-1));
		
		long [][] shiftsOut = new long[nImageN][];
		double [] listCC = new double[nImageN];
		
		double maxAverCC = (-1)*Double.MAX_VALUE;
		//final int dimShift = imageSet.imgs.get(0).numDimensions(); 
		int nIterMax = 0;
		ArrayList<long []> maxAverCCshifts = new ArrayList<long []>();
		
		for(i=0;i<nImageN;i++)
		{
			maxAverCCshifts.add(new long[nDimReg]);
		}
		
		//initial iteration
		int iterBeg = 0;
		//for zero iterations we want to do one without shifts
		if(nIterN == 0)
		{
			iterBeg=-1;
			limInterval = new FinalInterval(new long [nDimReg], new long [nDimReg]);
			normCC.lim_fractions = null;
			normCC.limInterval = limInterval;
			//no need to save, since they are the same
			bOutputInput = false;
		}
		
		long iterStartT, iterEndT;
		
		for(iter=iterBeg;(iter<nIterN && !bConverged);iter++)
		{
			IJ.showStatus("Averaging iteration "+Integer.toString(iter+1)+"...");
			avrgCC = 0.0;
			iterStartT = System.currentTimeMillis();
			//new average (sum and count)
			template.init(imgs_shift);
			//calculate shifts and CC values
			for(i=0;i<nImageN;i++)
			{
				
				//remove current image from the average
				currAverageImg = template.getTemplateForImage(imgs_shift.get(i));
				//currAverageImg = removeOneAverage(sumAndCount,imgs_shift.get(i));
				//ImageJFunctions.show(currAverageImg, "aver"+Integer.toString(i+1));
				//removeOneAverage
				normCC.caclulateGenNormCC(Views.zeroMin(currAverageImg), imageSet.imgs.get(i), false);
				
				avrgCC+=normCC.dMaxCC;
				listCC[i]=normCC.dMaxCC;
				if(nIterN >0)
				{
					IJ.showProgress(iter*(nImageN-1)+i,(nIterN)*(nImageN-1));
				}
				else
				{
					IJ.showProgress(i,nImageN);					
				}
				shiftsOut[i] = normCC.dShift.clone();
				for(j=0;j<normCC.dShift.length;j++)
				{
					shiftsOut[i][j]+=currAverageImg.min(j);
				}
				//shifts.set(i, normCC.dShift.clone());
	
			}
			// do median filtering of shifts
			medianCorrectShifts(shiftsOut);
			
			//see what we get and remember it
			avrgCC=avrgCC/nImageN;
			if(avrgCC > maxAverCC)
			{
				maxAverCC = avrgCC;
				nIterMax = iter+1;
				for(i=0;i<nImageN;i++)
				{
					for(d=0;d<nDimReg;d++)
					{
						maxAverCCshifts.get(i)[d] = shiftsOut[i][d];
					}
				}				
				
			}
			//put data to the table
			for(i=0;i<nImageN;i++)
			{
				shifts.set(i, shiftsOut[i]);
				//fill the table
				ptable.incrementCounter();
				ptable.addValue("iter", iter+1);
				ptable.addValue("particle", i+1);
				ptable.addValue("norm_CC_coeff", listCC[i]);
				for(k=0;k<normCC.dShift.length;k++)
				{
					ptable.addValue("coord_"+Integer.toString(k+1), shiftsOut[i][k]);
					ptable.addLabel(imageSet.image_names.get(i));
				}	
			}

			// calculate new img array with applied displacements			
			cumShift = buildShiftedIntervals(imageSet.imgs, imgs_shift, shifts);

			//sumAndCount = AverageWithoutZero.sumAndCountArray(imgs_shift);
			
			
			ptableCC.incrementCounter();
			ptableCC.addValue("iter", iter+1);
			ptableCC.addValue("averCC", avrgCC);
			ptableCC.addValue("cumShift", cumShift);
			iterEndT = System.currentTimeMillis();
			double elTime = (iterEndT -iterStartT)/60000.;
			String sTimeEl = " (time "; 
			if(elTime<1.0)
			{
				sTimeEl = sTimeEl +  df1.format(elTime*60.)+" sec)";
			}
			else
			{
				sTimeEl = sTimeEl +  df1.format(elTime)+" min)";
			}
			IJ.log("Iteration "+Integer.toString(iter+1)+" average CC " + df.format(avrgCC) +sTimeEl);			
			
			if(bSaveIntermediate && nIterN>0)
			{
				processIntermediate(iter+1);
			}
			
			//check if it is converged already
			if(Math.abs(oldAvrgCC - avrgCC)<0.000001 && Math.abs(oldCumShift - cumShift)<0.001)
			{
				bConverged = true;
			
				IJ.log("Converged before reaching the final iteration number");
			}
			else
			{
				oldCumShift = cumShift;
				oldAvrgCC = avrgCC;
			}

		}
		IJ.showStatus("Iterative averaging...done");
		
		IJ.showProgress(2,2);
		
		if(nIterN!=0)
		{
			IJ.log("Best result: iteration #" + Integer.toString(nIterMax) + " with average CC " + df.format(maxAverCC));
		
		}
		else
		{
			IJ.log("Iteration count is equal to zero, no registration was done, just averaging.");
		}
		ptable.show("Results");
		ptable.save(sPathGeneralOutput+"Results.csv");
		ptableCC.show("Average CC");
		ptableCC.save(sPathGeneralOutput+"Average_CC.csv");
		shifts = maxAverCCshifts;
		
		//calculate final average image		
		IJ.log("calculating final average image..");		
		imgs_avrg_out = getAlignedRAIs(shifts, true); 
		IntervalView<FloatType> finalAver = FinalOutput.averageArray(imgs_avrg_out, bIgnoreZeroInAveraging);
		String sOutFinalName = "";
		switch(nAveragingAim)
		{
			case TemplateAveraging.AVERAGE:
				sOutFinalName ="final_avg_";
			break;
			case TemplateAveraging.MASKED_AVERAGE:
				sOutFinalName ="final_mask_avg_";
			break;
		}
		sOutFinalName = sOutFinalName +Integer.toString(nIterMax);
		ImagePlus tempFin = MiscUtils.wrapFloatImgCal(finalAver,sOutFinalName,imageSet.cal,imageSet.bMultiCh);
		IJ.saveAsTiff(tempFin, sPathGeneralOutput+sOutFinalName);
		IJ.log("...done.");
		//calculate STD image
		if(bShowSD)
		{
			IJ.log("calculating final standard deviation image..");
			String sSDName = "final_std_"+Integer.toString(nIterMax);
			IntervalView<FloatType> finalSTD = FinalOutput.stdArray(imgs_avrg_out, finalAver, bIgnoreZeroInAveraging);
			ImagePlus finSD = MiscUtils.wrapFloatImgCal(finalSTD,sSDName,imageSet.cal,imageSet.bMultiCh);
			IJ.saveAsTiff(finSD, sPathGeneralOutput+sSDName);
			IJ.log("...done.");
		}
		if(bOutputInput)
		{

				IJ.log("saving registered inputs...");
				ImagePlus temp;
				for(i=0;i<nImageN;i++)
				{
					temp = MiscUtils.wrapFloatImgCal(Views.interval(Views.extendZero(imgs_avrg_out.get(i)),finalAver),"reg_"+imageSet.image_names.get(i), imageSet.cal, imageSet.bMultiCh); 
	
					IJ.saveAsTiff(temp, sPathRegistered+temp.getTitle());
				}
				IJ.log("..done");

			
		}
	}
	
	/** function calculating multi-channel shifts from imgs_multiCh given shifts 
	 * **/
	
	ArrayList<RandomAccessibleInterval< FloatType >> getMultiChAligned(final ArrayList<long []> shifts)//, String sTitle, Calibration cal)
	{
		final int nDim =  imageSet.imgs_multiCh.get(0).numDimensions();
		long [] curr_shift = new long [nDim];		
		
		final ArrayList<RandomAccessibleInterval< FloatType >> imgs_multiCh_reg = new ArrayList<RandomAccessibleInterval< FloatType >>(); 
		
		for(int iImCount=0; iImCount<shifts.size(); iImCount++)
		{
			//"jumping" over color channel, since it is xyczt and our shift is xyz			
			int j=0;
			for (int i=0;i<nDim;i++)
			{
				curr_shift[i]=shifts.get(iImCount)[j];
				if(i==1)
				{
					curr_shift[2]=0;
					i++;
				}
				j++;
			}			
			imgs_multiCh_reg.add(Views.translate(imageSet.imgs_multiCh.get(iImCount), curr_shift));
		}
		return imgs_multiCh_reg;
	}
	
	
	/** given input image array imgs_in and shifts, generates corresponding array of applied shifted interval views imgs_out **/
	public double buildShiftedIntervals(final ArrayList<RandomAccessibleInterval< FloatType >> imgs_in, final ArrayList<RandomAccessibleInterval< FloatType >> imgs_out, final ArrayList<long []> shifts)
	{
		// calculate new displacements
		double cumShift=0.0;
		double nCurrShift = 0.0;
		int i,j;
		int nDim = shifts.get(0).length;
		imgs_out.clear();
		for(i=0;i<imgs_in.size();i++)
		{
			imgs_out.add(i,Views.translate(imgs_in.get(i),shifts.get(i)));
			
			//estimate total displacement
			nCurrShift = 0.0;
			for(j=0;j<nDim;j++)
			{
				nCurrShift += shifts.get(i)[j]*shifts.get(i)[j];
			}
			cumShift+=Math.sqrt(nCurrShift);
		}
		return cumShift;
	}
	
	/** function generates shifts so all images are centered (+/- one pixel) if nMethod ==0 
	 * and just zero values if nMethod == 1 **/
	public ArrayList<long []> initShifts(final ArrayList<RandomAccessibleInterval< FloatType >> imgs_in, int nMethod)
	{
		int nDim = imgs_in.get(0).numDimensions();
		long [] currSize = imgs_in.get(0).dimensionsAsLongArray();
		long [] currMin = imgs_in.get(0).minAsLongArray();
		long [] currShift; 
		ArrayList<long []> shifts = new ArrayList<long []> ();		
		int i,j;
		int nDisp;
		for (i=0;i<imgs_in.size();i++)
		{
			
				imgs_in.get(i).dimensions(currSize);
				imgs_in.get(i).min(currMin);
				currShift = new long [nDim];
				if(nMethod == 0)
				{
					for(j=0;j<nDim;j++)
					{
						nDisp = -1*(int)Math.ceil(0.5*currSize[j]);
						currShift[j] = nDisp-currMin[j];
					}
				}
				shifts.add(currShift);

		}		
		return shifts;
	}
	
	public void medianCorrectShifts(final long [][] shifts_in)
	{
		int j;
		final int d = shifts_in[0].length;
		final int imN = shifts_in.length;
		final long [] dimX = new long[imN];
		long median;
		for(int i=0;i<d;i++)
		{
			for(j=0;j<imN;j++)
			{
				dimX[j]=shifts_in[j][i];
			}
			Arrays.sort(dimX);
			if (imN % 2 == 0)
			    median = Math.round(((double)dimX[imN/2] + (double)dimX[imN/2 - 1])/2.0);
			else
			    median = dimX[imN/2];
			for(j=0;j<imN;j++)
			{
				shifts_in[j][i]-=median;
			}
			
		}
		
	}
	
	public void processIntermediate(final int nIt)
	{
		ArrayList<RandomAccessibleInterval< FloatType >> imgs_avrg_out = getAlignedRAIs(shifts, true); 
		String sName =""; 
		
		switch(nAveragingAim)
		{
			case TemplateAveraging.AVERAGE:
				sName="it_avg_";
				break;
			case TemplateAveraging.MASKED_AVERAGE:
				sName="it_mask_avg_";
				break;
			//case TemplateAveraging.MEDIAN:
			//	sName="it_med_";
			//	break;
		}
		
		sName = sName+Integer.toString(nIt);
		ImagePlus temp;
		temp = MiscUtils.wrapFloatImgCal(FinalOutput.averageArray(imgs_avrg_out, bIgnoreZeroInAveraging),sName,imageSet.cal,imageSet.bMultiCh);
		if(bSaveIntermediate)
		{
			IJ.saveAsTiff(temp, sPathIntermediate +temp.getTitle());
		}
		else
		{
			temp.show();
		}	
	}
	
	public boolean dialogSettings()
	{
		int d;
		
		//double format formatting tool
		final DecimalFormatSymbols symbols = new DecimalFormatSymbols();
		symbols.setDecimalSeparator('.');
		final DecimalFormat df1 = new DecimalFormat ("#.#", symbols);
		
		sDims = imageSet.sRefDims;
		
		nDimReg = sDims.length();
		if(imageSet.bMultiCh)
		{
			nDimReg--; //remove the C component
		}
		
		
		limName = new Label[nDimReg];
		limVal = new TextField[nDimReg];
		double [] dLimits = new double [nDimReg];
		
		final String[] sIniTemplate = new String[ ]{"Centered","Zero (top-left)"};
		final String[] sAveragingAim = new String[ ]{"Average","Zero masked average"}; //"Median",
		final String[] limitsReg = new String[  ] {"No","by voxels", "by image fraction"};
		final GenericDialog gd1 = new GenericDialog( "Averaging parameters" );
		if(imageSet.bMultiCh)
		{
			final String[] channels = new String[ imageSet.nChannels];
			for ( int c = 0; c < channels.length; ++c )
				channels[ c ] = "use channel " + Integer.toString(c+1);
			gd1.addChoice( "For alignment", channels, channels[ 0 ] );
		}
		gd1.addChoice( "Initial template:", sIniTemplate, Prefs.get("RegisterNDFFT.IA.nIniTemplate", sIniTemplate[0]) );
		gd1.addNumericField("Number of iterations", Prefs.get("RegisterNDFFT.IA.nIterN",10),0);
		gd1.addChoice( "Template calculation:", sAveragingAim, Prefs.get("RegisterNDFFT.IA.nAveragingAim", sAveragingAim[0]) );
		gd1.addCheckbox("Use zero masked CC?", Prefs.get("RegisterNDFFT.IA.bExcludeZeros", false));		
		String sCurrChoice = Prefs.get("RegisterNDFFT.IA.sConstrain", "No");
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
					gd1.addNumericField(sDims.charAt(d)+" limit (px)", Prefs.get("RegisterNDFFT.IA.dMax"+sDims.charAt(d)+"px", 10.0), 3);
					break;
				case "by image fraction":
					gd1.addNumericField(sDims.charAt(d)+" limit (0-1)", Prefs.get("RegisterNDFFT.IA.dMax"+sDims.charAt(d)+"fr", 0.5), 3);
					break;
					
			}
			limName[d] = gd1.getLabel();
			limVal[d] = (TextField)gd1.getNumericFields().get(d+1);	
			if(sCurrChoice.equals("No"))
			{
				limVal[d].setEnabled(false);
			}
		}
		gd1.addMessage("Output:");
		gd1.addCheckbox("Save intermediates?", Prefs.get("RegisterNDFFT.IA.bSaveIntermediate",false));
		gd1.addCheckbox("Save registered inputs?", Prefs.get("RegisterNDFFT.IA.bOutputInput",false));
		gd1.addCheckbox("Save SD image?", Prefs.get("RegisterNDFFT.IA.bShowSD",false));

		gd1.addDialogListener(this);
		gd1.showDialog();
		
		if ( gd1.wasCanceled() )
			return false;
		
		if(imageSet.bMultiCh)
		{
			imageSet.alignChannel = gd1.getNextChoiceIndex();
		}
		
		nIniTemplate = gd1.getNextChoiceIndex();
		Prefs.set("RegisterNDFFT.IA.nIniTemplate", sIniTemplate[nIniTemplate]);
		
		nIterN  = (int)gd1.getNextNumber();
		Prefs.set("RegisterNDFFT.IA.nIterN", nIterN);
		
		nAveragingAim = gd1.getNextChoiceIndex();
		Prefs.set("RegisterNDFFT.IA.nAveragingAim", sAveragingAim[nAveragingAim]);
		
		bZeroMask  = gd1.getNextBoolean();
		Prefs.set("RegisterNDFFT.IA.bExcludeZeros", bZeroMask);
		
		nConstrainReg = gd1.getNextChoiceIndex();
		Prefs.set("RegisterNDFFT.IA.sConstrain", limitsReg[nConstrainReg]);
		
		IJ.log("Initial template: "+ sIniTemplate[nIniTemplate]);
		
		IJ.log("Template calculation: "+ sAveragingAim[nAveragingAim]);
		if(nAveragingAim == TemplateAveraging.AVERAGE)
		{
			bIgnoreZeroInAveraging = false;
		}
		if(nAveragingAim == TemplateAveraging.MASKED_AVERAGE)
		{
			bIgnoreZeroInAveraging = true;
		}

		
		if(bZeroMask)
		{
			IJ.log("Using zero masked cross-correlation.");
		}
		else
		{
			IJ.log("Using non-masked cross-correlation.");			
		}
		
		
		if(nIterN==0)
		{
			IJ.log("Iteration count is equal to zero, no registration, only CC/average calculation.");
		}
		else
		{
		
			if(nConstrainReg == 0)
			{
				IJ.log("Averaging without constrains.");
			}
			else
			{
				if(nConstrainReg == 1)
				{
					IJ.log("Averaging with constrain specified in voxels:");
	
					for(d=0;d<nDimReg;d++)
					{
						dLimits[d]=Math.abs(gd1.getNextNumber());
						Prefs.set("RegisterNDFFT.IA.dMax"+sDims.charAt(d)+"px",dLimits[d]);
						IJ.log("Axis " +sDims.charAt(d)+": "+df1.format(dLimits[d])+" pixels");
					}
					
				}
				else
				{
					IJ.log("Averaging with constrain specified as a fraction of max displacement:");
					for(d=0;d<nDimReg;d++)
					{
						dLimits[d]=Math.min(Math.abs(gd1.getNextNumber()), 1.0);
						Prefs.set("RegisterNDFFT.IA.dMax"+sDims.charAt(d)+"fr",dLimits[d]);
						IJ.log("Axis " +sDims.charAt(d)+": "+df1.format(dLimits[d]));
					}
				}
			}
			
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
		}
		
		bSaveIntermediate = gd1.getNextBoolean();
		Prefs.set("RegisterNDFFT.IA.bSaveIntermediate", bSaveIntermediate);
		
		bOutputInput = gd1.getNextBoolean();
		Prefs.set("RegisterNDFFT.IA.bOutputInput", bOutputInput);
		
		bShowSD = gd1.getNextBoolean();
		Prefs.set("RegisterNDFFT.IA.bShowSD", bShowSD);
		
		DirectoryChooser dc = new DirectoryChooser ( "Choose a folder to save output..." );
		sPathGeneralOutput = dc.getDirectory();
		if(sPathGeneralOutput == null)
			return false;
		
		if(bSaveIntermediate && nIterN>0)
		{
			sPathIntermediate  = sPathGeneralOutput +"intermediate/";
			File theDir = new File(sPathIntermediate);
			if (!theDir.exists())
			{
			    theDir.mkdirs();
			}
			
		}
		if(bOutputInput)
		{
			sPathRegistered  = sPathGeneralOutput +"registered/";
			File theDir = new File(sPathRegistered);
			if (!theDir.exists())
			{
			    theDir.mkdirs();
			}
		}
		return true;

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
							limVal[d].setText(df1.format(Prefs.get("RegisterNDFFT.IA.dMax"+sDims.charAt(d)+"px", 10.0)));
						}
						break;
					case 2:
						for(d=0;d<nDimReg;d++)
						{
							limName[d].setText(sDims.charAt(d)+" limit (0-1)");
							limVal[d].setEnabled(true);
							limVal[d].setText(df1.format(Prefs.get("RegisterNDFFT.IA.dMax"+sDims.charAt(d)+"fr", 0.5)));

						}
						break;
				}
			}
		}
		return true;
	}
	
	/** function returns aligned intervals **/
	ArrayList<RandomAccessibleInterval< FloatType >> getAlignedRAIs(final ArrayList<long []> shifts_, boolean bIgnoreZero)
	{
		ArrayList<RandomAccessibleInterval< FloatType >> imgs_avrg_out; 
		
		if(imageSet.bMultiCh)
		{
			imgs_avrg_out = getMultiChAligned(shifts_);
		}
		else
		{
			// calculate new img array with applied displacements	
			imgs_avrg_out = new ArrayList<RandomAccessibleInterval< FloatType >>();
			buildShiftedIntervals(imageSet.imgs, imgs_avrg_out,shifts_);
		}
		return imgs_avrg_out;
	}
		
	public static void main( final String[] args )
	{
		// open an ImageJ window
		 new ImageJ();
		//IJ.open("/home/eugene/Desktop/projects/RegisterNDFFT/single/MAX_089-1.tif");
		//IJ.open("/home/eugene/Desktop/projects/RegisterNDFFT/single/MAX_098-1.tif");
		IterativeAveraging it = new IterativeAveraging();
		it.run(null);

	}


}
