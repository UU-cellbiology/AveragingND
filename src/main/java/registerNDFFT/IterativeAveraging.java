package registerNDFFT;

import java.io.IOException;

import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.io.DirectoryChooser;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import io.scif.img.ImgIOException;
import io.scif.img.ImgOpener;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.exception.IncompatibleTypeException;
import net.imglib2.img.ImagePlusAdapter;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;

public class IterativeAveraging implements PlugIn {

	public int nIterN = 0;
	public int nImageN = 0;
	public double dMaxFraction = 0.5;
	public int alignChannel = 1;
	public boolean bShowIntermediateAverage=false;
	public boolean bMultiCh = false;
	public int nIniTemplate = 0;
	public int nInput = 0;
	public boolean bExcludeZeros = false;
	public boolean bOutputInput = false;
	
	public int numChannels = 1;
	
	/** original images (with full channels) **/
	ArrayList<RandomAccessibleInterval< FloatType >> imgs_multiCh;
	
	/** only channel used for alignment **/
	ArrayList<RandomAccessibleInterval< FloatType >> imgs;
	
	/** shifted images for average calculation**/
	ArrayList<RandomAccessibleInterval< FloatType >> imgs_shift;
	
	/** shift of each image **/
	ArrayList<long []> shifts;
	
	/** shift of each image **/
	ArrayList<String> image_names;
	
	Calibration calibInput;
	
	@Override
	public void run(String paramString) {
		// TODO Auto-generated method stub
		int i,j,k, iter;
	
		
		//final String[] sIniTemplate = new String[3];
		final String[] sIniTemplate = new String[2];
		sIniTemplate[0] = "Average (center)";
		sIniTemplate[1] = "Average (top-left)";
		//sIniTemplate[2] = "Least squares (pairwise)";
		
		final String[] sInput = new String[2];
		sInput[0] = "All currently open images";
		sInput[1] = "Specify images in a folder";
		
		final GenericDialog gd = new GenericDialog( "Iterative registration" );

		gd.addChoice( "Input images:", sInput, Prefs.get("RegisterNDFFT.IA.nInput", sInput[0]) );
		gd.addChoice( "Initial template:", sIniTemplate, Prefs.get("RegisterNDFFT.IA.nIniTemplate", sIniTemplate[0]) );
		gd.addNumericField("Number of iterations", Prefs.get("RegisterNDFFT.IA.nIterN",10),0);
		gd.addCheckbox("Exclude zero values?", Prefs.get("RegisterNDFFT.IA.bExcludeZeros", false));	
		gd.addNumericField("Maximum shift (fraction, 0-1 range)", Prefs.get("RegisterNDFFT.IA.dMaxFraction", 0.4), 3);
		gd.addCheckbox("Show intermediate average", Prefs.get("RegisterNDFFT.IA.bShowIntermediateAverage",false));
		gd.addCheckbox("Output registered inputs?", Prefs.get("RegisterNDFFT.IA.bOutputInput",false));
		gd.showDialog();
		
		if ( gd.wasCanceled() )
			return;				

		nInput = gd.getNextChoiceIndex();
		Prefs.set("RegisterNDFFT.IA.nInput", sInput[nInput]);
		nIniTemplate = gd.getNextChoiceIndex();
		Prefs.set("RegisterNDFFT.IA.nIniTemplate", sIniTemplate[nIniTemplate]);
		nIterN  = (int)gd.getNextNumber();
		Prefs.set("RegisterNDFFT.IA.nIterN", nIterN);
		bExcludeZeros  = gd.getNextBoolean();
		Prefs.set("RegisterNDFFT.IA.bExcludeZeros", bExcludeZeros);
		dMaxFraction  = gd.getNextNumber();
		Prefs.set("RegisterNDFFT.IA.dMaxFraction", dMaxFraction);
		bShowIntermediateAverage = gd.getNextBoolean();
		Prefs.set("RegisterNDFFT.IA.bShowIntermediateAverage", bShowIntermediateAverage);
		bOutputInput = gd.getNextBoolean();
		Prefs.set("RegisterNDFFT.IA.bOutputInput", bOutputInput);
		
		//init image arrays		
		imgs_multiCh = new ArrayList<RandomAccessibleInterval< FloatType >>();
		imgs = new ArrayList<RandomAccessibleInterval< FloatType >>();		
		imgs_shift = new ArrayList<RandomAccessibleInterval< FloatType >>();
		shifts = new ArrayList<long []>();
		image_names = new ArrayList<String>();
		if(nInput == 0)
		{
			 if(!loadAllOpenImages())
				 return;			 
		}
		else
		{
			if(!loadFolderTiff())
				return;
		}
		
		// case nIniTemplate ==1 automatically happens during image loading
		//for ==0 we need to generate new array
		if(nIniTemplate == 0)
		{
			shifts = centeredShifts(imgs);
		}
		//create a new shifted array of images
		buildShiftedIntervals(imgs, imgs_shift,shifts);
		

		ArrayList<IntervalView<FloatType>> sumAndCount = null;
		sumAndCount = AverageWithoutZero.sumAndCountArray(imgs_shift);

		IntervalView<FloatType> currAverageImg;
		if(bShowIntermediateAverage)
		{		
			MiscUtils.wrapFloatImgCal(AverageWithoutZero.averageFromSumAndCount(sumAndCount), "average iteration 0", calibInput, false).show();
		}
		
		GenNormCC normCC = new GenNormCC();
		normCC.bVerbose = false;
		normCC.bExcludeZeros=bExcludeZeros;

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
		IJ.showProgress(0, (nIterN)*(nImageN-1));
		
		long [][] shiftsOut = new long[nImageN][];
		double [] listCC = new double[nImageN];
		
		double maxAverCC = (-1)*Double.MAX_VALUE;
		int nIterMax = 0;
		ArrayList<long []> maxAverCCshifts = new ArrayList<long []>();
		final int dimShift= imgs.get(0).numDimensions();
		
		for(i=0;i<nImageN;i++)
		{
			maxAverCCshifts.add(new long[dimShift]);
		}
		
		DecimalFormatSymbols symbols = new DecimalFormatSymbols();
		symbols.setDecimalSeparator('.');
		DecimalFormat df = new DecimalFormat ("#.########", symbols);
		DecimalFormat df1 = new DecimalFormat ("#.#", symbols);
		
		long iterStartT, iterEndT;
		
		for(iter=0;(iter<nIterN && !bConverged);iter++)
		{
			IJ.showStatus("Averaging iteration "+Integer.toString(iter+1)+"...");
			avrgCC = 0.0;
			iterStartT = System.currentTimeMillis();
			//calculate shifts and CC values
			for(i=0;i<nImageN;i++)
			{
				//remove current image from the average
				currAverageImg = removeOneAverage(sumAndCount,imgs_shift.get(i));
				//ImageJFunctions.show(currAverageImg, "aver"+Integer.toString(i+1));
				//removeOneAverage
				normCC.caclulateGenNormCC(Views.zeroMin(currAverageImg), imgs.get(i), dMaxFraction , false);
				
				avrgCC+=normCC.dMaxCC;
				listCC[i]=normCC.dMaxCC;
				IJ.showProgress(iter*(nImageN-1)+i,(nIterN)*(nImageN-1));
				
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
					for(int d=0;d< dimShift;d++)
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
					ptable.addLabel(image_names.get(i));
				}	
			}

			// calculate new img array with applied displacements			
			cumShift = buildShiftedIntervals(imgs, imgs_shift,shifts);
			//new average (sum and count)
			sumAndCount = AverageWithoutZero.sumAndCountArray(imgs_shift);
			
			
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
			
			if(bShowIntermediateAverage)
			{
				MiscUtils.wrapFloatImgCal(AverageWithoutZero.averageFromSumAndCount(sumAndCount),"average iteration "+Integer.toString(iter+1),calibInput, false).show();
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
		
			ptable.show("Results");
			ptableCC.show("Average CC");
			shifts = maxAverCCshifts;
		}
		else
		{
			IJ.log("Iteration count is equal to zero, no registration was done, just averaging.");
		}
		

		
		ArrayList<RandomAccessibleInterval< FloatType >> imgs_avrg_out = new ArrayList<RandomAccessibleInterval< FloatType >>();
		
		
		
		if(bMultiCh)
		{
			getMultiChAligned(imgs_avrg_out, shifts);
		}
		else
		{
			// calculate new img array with applied displacements	
			buildShiftedIntervals(imgs, imgs_avrg_out,shifts);
			// = imgs_shift;
		}
		
		IJ.log("generating average image..");
		
		//calculate final average image
		IntervalView<FloatType> finalAver = AverageWithoutZero.averageArray(imgs_avrg_out, true);
		IJ.log("...done.");
		MiscUtils.wrapFloatImgCal(finalAver,"final_average_"+Integer.toString(nIterMax),calibInput,bMultiCh).show();
		
		//calculate STD image
		IJ.log("generating standard deviation image..");
		IntervalView<FloatType> finalSTD = AverageWithoutZero.stdArray(imgs_avrg_out, finalAver, true);
		MiscUtils.wrapFloatImgCal(finalSTD,"final_std_"+Integer.toString(nIterMax),calibInput,bMultiCh).show();
		IJ.log("...done.");
		
		if(bOutputInput)
		{
				DirectoryChooser dc = new DirectoryChooser ( "Choose a folder to save output..." );
				String sPath = dc.getDirectory();
				if(sPath == null || finalAver == null)
					return;
				ImagePlus temp;
				for(i=0;i<nImageN;i++)
				{
					temp = MiscUtils.wrapFloatImgCal(Views.interval(Views.extendZero(imgs_avrg_out.get(i)),finalAver),"iter_aver_"+image_names.get(i),calibInput, bMultiCh); 
	
					IJ.saveAsTiff(temp, sPath+temp.getTitle());
				}
			
		}
	}
	
	/** function calculating multi-channel shifts from imgs_multiCh given shifts 
	 * **/
	//void showMultiChAverage(final ArrayList<RandomAccessibleInterval< FloatType >> imgs_multiCh, final ArrayList<long []> shifts, String sTitle, Calibration cal)
	void getMultiChAligned(final ArrayList<RandomAccessibleInterval< FloatType >> imgs_multiCh_reg, final ArrayList<long []> shifts)//, String sTitle, Calibration cal)
	{
		int nDim = imgs_multiCh.get(0).numDimensions();
		long [] curr_shift = new long [nDim];
		int iImCount;
		int i,j;
		//final ArrayList<RandomAccessibleInterval< FloatType >> imgs_multiCh_reg =new ArrayList<RandomAccessibleInterval< FloatType >>(); 
		
		for(iImCount=0;iImCount<shifts.size();iImCount++)
		{

			//"jumping" over color channel, since it is xyczt and our shift is xyz
			
			j=0;
			for (i=0;i<nDim;i++)
			{
				curr_shift[i]=shifts.get(iImCount)[j];
				if(i==1)
				{
					curr_shift[2]=0;
					i++;
				}
				j++;
			}
			
			imgs_multiCh_reg.add(Views.translate(imgs_multiCh.get(iImCount), curr_shift));
		}
		
	}
	
	/** given Sum and Count images alSumCnt, this function subtracts removedImage from Sum,
	 * reduces corresponding Count voxels and returns averaged image (with coordinate origin at (0, 0, ..0) **/	
	IntervalView<FloatType> removeOneAverage(ArrayList<IntervalView<FloatType>> alSumCnt, RandomAccessibleInterval< FloatType > removedImage)
	{
	
		long [] origin = alSumCnt.get(0).minAsLongArray();
		final Img<FloatType> avrgImgArr = ArrayImgs.floats(alSumCnt.get(0).dimensionsAsLongArray());
		final IntervalView<FloatType> avrgImg = Views.translate(avrgImgArr, origin );		
		final IntervalView<FloatType> removeInt = Views.interval(Views.extendZero(removedImage), avrgImg);
		Cursor<FloatType> avrgC = avrgImg.cursor();
		Cursor<FloatType> remC = removeInt.cursor();
		Cursor<FloatType> sumC = alSumCnt.get(0).cursor();
		Cursor<FloatType> cntC = alSumCnt.get(1).cursor();
		float fCnt, fRem, fSum;
		while(avrgC.hasNext())
		{
			
			avrgC.fwd();
			remC.fwd();
			sumC.fwd();
			cntC.fwd();
			fSum = sumC.get().get();
			fRem = remC.get().get();
			fCnt = cntC.get().get();
			if(fRem>0.0f)
			{
				fSum-=fRem;
				fCnt--;
			}		
			
			if (fCnt>0.5f)
			{
				avrgC.get().set(fSum/fCnt);
			}
			else
			{
				avrgC.get().set(0.0f);
			}
		}
				
		//return Views.zeroMin(avrgImg);
		return avrgImg;
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
	/** function generates shifts so all images are centered (+/- one pixel) **/
	public ArrayList<long []> centeredShifts(final ArrayList<RandomAccessibleInterval< FloatType >> imgs_in)
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
			for(j=0;j<nDim;j++)
			{
				nDisp = -1*(int)Math.ceil(0.5*currSize[j]);
				currShift[j] = nDisp-currMin[j];
			}
			shifts.add(currShift);
		}		
		return shifts;
	}
	/** function fills analysis arrays with images currently open in ImageJ **/
	public boolean loadAllOpenImages()
	{
		int i;
		final int[] idList = WindowManager.getIDList();		

		if ( idList == null || idList.length < 2 )
		{
			IJ.error( "You need at least two open images." );
			return false;
		}
		
		
		calibInput = WindowManager.getImage(idList[0]).getCalibration();		
		nImageN=idList.length;
		
		numChannels = WindowManager.getImage( idList[0] ).getNChannels();
		
		if(numChannels>1)
		{
			bMultiCh = true;
			final GenericDialog gdCh = new GenericDialog( "Choose registration channel" );
			
			final String[] channels = new String[ numChannels ];
			for ( int c = 0; c < channels.length; ++c )
				channels[ c ] = "use channel " + Integer.toString(c+1);
			gdCh.addChoice( "For alignment", channels, channels[ 0 ] );
			gdCh.showDialog();
			
			if ( gdCh.wasCanceled() )
				return false;				

			alignChannel = gdCh.getNextChoiceIndex();
		}

		for(i=0;i<nImageN;i++)
		{
			if(!bMultiCh)
			{
				imgs.add(ImagePlusAdapter.convertFloat(WindowManager.getImage(idList[i])));
			}
			else
			{
				imgs_multiCh.add(ImagePlusAdapter.convertFloat(WindowManager.getImage(idList[i])));
				imgs.add(Views.hyperSlice(imgs_multiCh.get(i),2,alignChannel));
			}
			shifts.add(new long [imgs.get(0).numDimensions()]);
			image_names.add(WindowManager.getImage(idList[i]).getTitle());
		}
		if(!bMultiCh)
		{
			IJ.log("Averaging "+ Integer.toString(nImageN) + " images.");
		}
		else
		{
			IJ.log("Averaging "+Integer.toString(nImageN) + " images with " + Integer.toString(numChannels) + " channels.");
			IJ.log("Using channel "+Integer.toString(alignChannel+1) + " for alignment.");
		}
		return true;
	}
	
	/** function opens folder with Tiff images and loads them to analysis arrays **/
	public boolean loadFolderTiff()
	{
		
		int i;
		DirectoryChooser dc = new DirectoryChooser ( "Choose a folder with images.." );
		String sPath = dc.getDirectory();
		List<String> files;
		if (sPath != null)
		{
			
			//IJ.log(sPath);
			try {

				files = MiscUtils.findFiles(Paths.get(sPath), ".tif");
				//files.forEach(x -> IJ.log(x));
			} catch (IOException e) {
				e.printStackTrace();
				return false;
            }
			nImageN=files.size();
			if (nImageN == 0)
				return false;
			ImagePlus impBridge =IJ.openImage(files.get(0));
			numChannels = impBridge.getNChannels();
			calibInput =  impBridge.getCalibration();
			
			if(numChannels>1)
			{
				bMultiCh = true;
				final GenericDialog gdCh = new GenericDialog( "Choose registration channel" );
				
				final String[] channels = new String[ numChannels ];
				for ( int c = 0; c < channels.length; ++c )
					channels[ c ] = "use channel " + Integer.toString(c+1);
				gdCh.addChoice( "For alignment", channels, channels[ 0 ] );
				gdCh.showDialog();
				
				if ( gdCh.wasCanceled() )
					return false;				

				alignChannel = gdCh.getNextChoiceIndex();
			}
			
			IJ.showProgress(0, nImageN-1);
			IJ.showStatus("Loading images..");		
			for(i=0;i<nImageN;i++)
			{
				impBridge = IJ.openImage(files.get(i));
			
				if(!bMultiCh)
				{
					imgs.add(ImagePlusAdapter.convertFloat(impBridge));
				}
				else
				{
					imgs_multiCh.add(ImagePlusAdapter.convertFloat(impBridge));
					imgs.add(Views.hyperSlice(imgs_multiCh.get(i),2,alignChannel));
				}
				shifts.add(new long [imgs.get(0).numDimensions()]);
				image_names.add(impBridge.getTitle());
				IJ.showProgress(i, nImageN-1);
			}
			IJ.showProgress(2,2);
			IJ.showStatus("Loading images..done.");	
			if(!bMultiCh)
			{
				IJ.log("Averaging "+ Integer.toString(nImageN) + " images.");
			}
			else
			{
				IJ.log("Averaging "+Integer.toString(nImageN) + " images with " + Integer.toString(numChannels) + " channels.");
				IJ.log("Using channel "+Integer.toString(alignChannel+1) + " for alignment.");
			}
		}
		
		else
		{
			return false;
		}
	
		
		return true;
	}
	
	public void medianCorrectShifts(final long [][] shifts_in)
	{
		int i,j;
		int d = shifts_in[0].length;
		int imN = shifts_in.length;
		long [] dimX = new long[imN];
		long median;
		for(i=0;i<d;i++)
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
	
	public static void main( final String[] args ) throws ImgIOException, IncompatibleTypeException
	{
		// open an ImageJ window
		new ImageJ();
		IterativeAveraging it = new IterativeAveraging();
		it.run(null);

		
	
	}
	
}
