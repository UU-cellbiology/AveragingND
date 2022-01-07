package fiji.plugin.RegisterNDFFT;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.io.DirectoryChooser;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import io.scif.config.SCIFIOConfig;
import io.scif.config.SCIFIOConfig.ImgMode;
import io.scif.img.ImgOpener;
import net.imglib2.Cursor;
import net.imglib2.RandomAccessibleInterval;
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
	
	public int numChannels = 1;
	
	/** original images (with full channels) **/
	ArrayList<RandomAccessibleInterval< FloatType >> imgs_multiCh;
	
	/** only channel used for alignment **/
	ArrayList<RandomAccessibleInterval< FloatType >> imgs;
	
	/** shifted images for average calculation**/
	ArrayList<RandomAccessibleInterval< FloatType >> imgs_shift;
	
	/** shift of each image **/
	ArrayList<long []> shifts;
	
	Calibration calibInput;
	
	@Override
	public void run(String paramString) {
		// TODO Auto-generated method stub
		int i,k, iter;
	
		
		final String[] sIniTemplate = new String[3];
		sIniTemplate[0] = "Average (center)";
		sIniTemplate[1] = "Average (top-left)";
		sIniTemplate[2] = "Least squares (pairwise)";
		
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

		
		//init image arrays		
		imgs_multiCh = new ArrayList<RandomAccessibleInterval< FloatType >>();
		imgs = new ArrayList<RandomAccessibleInterval< FloatType >>();		
		imgs_shift = new ArrayList<RandomAccessibleInterval< FloatType >>();
		shifts = new ArrayList<long []>();
		
		if(nInput ==0)
		{
			 if(!loadAllOpenImages())
				 return;			 
		}
		else
		{
			if(!loadFolderTiff())
				return;
		}
		
		
		//determine initial template
		if(nIniTemplate==0)
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
		double oldAvrgCC = -100.0;
		
		boolean bConverged = false;
		
		IJ.showProgress(0, (nIterN)*(nImageN-1));
		for(iter=0;(iter<nIterN && !bConverged);iter++)
		{
			IJ.showStatus("Averaging iteration "+Integer.toString(iter)+"...");
			avrgCC=0.0;
			
			for(i=0;i<nImageN;i++)
			{
				//remove current image from the average
				currAverageImg = removeOneAverage(sumAndCount,imgs_shift.get(i));

				normCC.caclulateGenNormCC(currAverageImg, imgs.get(i), dMaxFraction , false);
				ptable.incrementCounter();
				ptable.addValue("iter", iter+1);
				ptable.addValue("particle", i+1);
				ptable.addValue("norm_CC_coeff", normCC.dMaxCC);
				for(k=0;k<normCC.dShift.length;k++)
				{
					ptable.addValue("coord_"+Integer.toString(k+1), normCC.dShift[k]);
				}
				
				shifts.set(i, normCC.dShift.clone());
				avrgCC+=normCC.dMaxCC;
				IJ.showProgress(iter*(nImageN-1)+i,(nIterN)*(nImageN-1));
			}
			
			// calculate new img array with applied displacements			
			cumShift= buildShiftedIntervals(imgs, imgs_shift,shifts);
			//new average (sum and count)
			sumAndCount = AverageWithoutZero.sumAndCountArray(imgs_shift);
			
			avrgCC=avrgCC/nImageN;
			ptableCC.incrementCounter();
			ptableCC.addValue("iter", iter+1);
			ptableCC.addValue("averCC", avrgCC);
			ptableCC.addValue("cumShift", cumShift);
			IJ.log("Iteration "+Integer.toString(iter+1)+" average CC " +Double.toString(avrgCC));
			
			//averageImg = AverageWithoutZero.averageArray(imgs_shift);
			if(bShowIntermediateAverage)
			{
				MiscUtils.wrapFloatImgCal(AverageWithoutZero.averageFromSumAndCount(sumAndCount),"average iteration "+Integer.toString(iter+1),calibInput, false).show();
			}
			
			//check if it is converged already
			if(Math.abs(oldAvrgCC - avrgCC)<0.000001 && Math.abs(oldCumShift - cumShift)<0.001)
			{
				bConverged = true;
				iter--;
			}
			else
			{
				oldCumShift = cumShift;
				oldAvrgCC = avrgCC;
			}

		}
		IJ.showProgress((nIterN-1)*(nImageN),(nIterN)*(nImageN-1));
		ptable.show("Results");
		ptableCC.show("Average CC");
		if(!bMultiCh)
		{
			MiscUtils.wrapFloatImgCal(AverageWithoutZero.averageFromSumAndCount(sumAndCount),"final average "+Integer.toString(iter+1),calibInput,false).show();
		}
		else
		{
			showMultiChAverage(imgs_multiCh,shifts,"final average "+Integer.toString(iter+1),calibInput);
		}
	}
	
	/** function calculating multi-channel average of imgs_multiCh given shifts 
	 * and diplays it as a proper ImagePlus **/
	void showMultiChAverage(final ArrayList<RandomAccessibleInterval< FloatType >> imgs_multiCh, final ArrayList<long []> shifts, String sTitle, Calibration cal)
	{
		int nDim = imgs_multiCh.get(0).numDimensions();
		long [] curr_shift = new long [nDim];
		int iImCount;
		int i,j;
		ArrayList<RandomAccessibleInterval< FloatType >> imgs_multiCh_reg =new ArrayList<RandomAccessibleInterval< FloatType >>(); 
		
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
		IntervalView<FloatType>  averageImg = AverageWithoutZero.averageArray(imgs_multiCh_reg);
		MiscUtils.wrapFloatImgCal(averageImg,sTitle,cal,true).show();
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
				
		return Views.zeroMin(avrgImg);
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

				files = MiscUtils.findFiles(Paths.get(sPath), "tif");
				//files.forEach(x -> IJ.log(x));
			} catch (IOException e) {
				e.printStackTrace();
				return false;
            }
			nImageN=files.size();
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
	

}
