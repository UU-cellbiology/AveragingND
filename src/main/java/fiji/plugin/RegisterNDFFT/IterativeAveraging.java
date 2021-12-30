package fiji.plugin.RegisterNDFFT;

import java.util.ArrayList;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.ImagePlusAdapter;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;

public class IterativeAveraging implements PlugIn {

	public int nIterN;
	public int nImageN;
	public double dMaxFraction = 0.5;
	public int alignChannel = 1;
	public boolean bShowIntermediateAverage=false;
	public boolean bMultiCh = false;
	
	@Override
	public void run(String paramString) {
		// TODO Auto-generated method stub
		int i,j,k, iter;
		int nDim;
		//this part is honestly stolen from "Pairwise Stitching" plugin
		//https://github.com/fiji/Stitching/blob/master/src/main/java/plugin/Stitching_Pairwise.java
		// get list of image stacks
		final int[] idList = WindowManager.getIDList();		

		if ( idList == null || idList.length < 2 )
		{
			IJ.error( "You need at least two open images." );
			return;
		}
		
		// create channel selector		
		ImagePlus imp = WindowManager.getImage( idList[0] );		
		final int numChannels = imp.getNChannels();
		final String[] channels = new String[ numChannels ];
		
		final GenericDialog gd = new GenericDialog( "Iterative registration" );
		if(numChannels>1)
		{
			bMultiCh = true;
			for ( int c = 0; c < channels.length; ++c )
				channels[ c ] = "use channel " + Integer.toString(c+1);
			gd.addChoice( "For alignment", channels, channels[ 0 ] );
		}

		gd.addNumericField("Number of iterations", Prefs.get("RegisterNDFFT.IA.nIterN",10),0);
		gd.addNumericField("Maximum shift (fraction, 0-1 range)", Prefs.get("RegisterNDFFT.IA.dMaxFraction", 0.4), 3);
		gd.addCheckbox("Show intermediate average", Prefs.get("RegisterNDFFT.IA.bShowIntermediateAverage",false));
		gd.showDialog();
		
		if ( gd.wasCanceled() )
			return;
		
		if(bMultiCh)
		{
			alignChannel = gd.getNextChoiceIndex();
		}
		nIterN  = (int)gd.getNextNumber();
		Prefs.set("RegisterNDFFT.IA.nIterN", nIterN);
		dMaxFraction  = gd.getNextNumber();
		Prefs.set("RegisterNDFFT.IA.dMaxFraction", dMaxFraction);
		bShowIntermediateAverage = gd.getNextBoolean();
		Prefs.set("RegisterNDFFT.IA.bShowIntermediateAverage", bShowIntermediateAverage);

		
		ArrayList<RandomAccessibleInterval< FloatType >> imgs_multiCh = new ArrayList<RandomAccessibleInterval< FloatType >>();
		/** original images (or channels of images) **/
		ArrayList<RandomAccessibleInterval< FloatType >> imgs = new ArrayList<RandomAccessibleInterval< FloatType >>();
		/** shifted images for average calculation**/
		ArrayList<RandomAccessibleInterval< FloatType >> imgs_shift = new ArrayList<RandomAccessibleInterval< FloatType >>();
		/** shift of each image **/
		ArrayList<long []> shifts = new ArrayList<long[]>();
		nImageN=idList.length;
		
		Calibration cal;
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
			//imgs_shift.add(Views.translate(imgs.get(i),shifts.get(i)));
			imgs_shift.add(imgs.get(i));
		}
		nDim = imgs.get(0).numDimensions();
		cal=WindowManager.getImage(idList[0]).getCalibration();
		
		ArrayList<IntervalView<FloatType>> sumAndCount = null;
		sumAndCount = AverageWithoutZero.sumAndCountArray(imgs);
		//IntervalView<FloatType> averageImg = AverageWithoutZero.averageArray(imgs);
		IntervalView<FloatType> averageImg = null;// AverageWithoutZero.averageArray(imgs);
		Img<FloatType> currAverageImg;// = AverageWithoutZero.averageArray(imgs);
		if(bShowIntermediateAverage)
		{
			//MiscUtils.wrapFloatImgCal(averageImg, "average 0",cal, false).show();
			
			MiscUtils.wrapFloatImgCal(AverageWithoutZero.averageFromSumAndCount(sumAndCount), "average 0",cal, false).show();
		}
		
		GenNormCC normCC = new GenNormCC();
		double avrgCC;
		ResultsTable ptable = ResultsTable.getResultsTable();
		ResultsTable ptableCC = new ResultsTable();
		ptableCC.setPrecision(8);
		double cumShift=0.0;
		double nCurrShift=0.0;
		IJ.showProgress(0, (nIterN)*(nImageN-1));
		for(iter=0;iter<nIterN;iter++)
		{

			avrgCC=0.0;
			
			for(i=0;i<nImageN;i++)
			{
				//TODO: proper subtraction
				currAverageImg = subtractFromAverage(averageImg,imgs_shift.get(i));
				ImageJFunctions.show(currAverageImg).setTitle( "sep" );
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
			cumShift=0.0;
			
			for(i=0;i<nImageN;i++)
			{
				imgs_shift.set(i,Views.translate(imgs.get(i),shifts.get(i)));
				
				//estimate total displacement
				nCurrShift = 0.0;
				for(j=0;j<nDim;j++)
				{
					nCurrShift += shifts.get(i)[j]*shifts.get(i)[j];
				}
				cumShift+=Math.sqrt(nCurrShift);
			}
			avrgCC=avrgCC/nImageN;
			ptableCC.incrementCounter();
			ptableCC.addValue("iter", iter+1);
			ptableCC.addValue("averCC", avrgCC);
			ptableCC.addValue("cumShift", cumShift);
			IJ.log("Iteration "+Integer.toString(iter+1)+" average CC " +Double.toString(avrgCC));
			averageImg = AverageWithoutZero.averageArray(imgs_shift);
			if(bShowIntermediateAverage)
			{
				MiscUtils.wrapFloatImgCal(averageImg,"average iteration "+Integer.toString(iter+1),cal, false).show();
			}

		}
		IJ.showProgress((nIterN-1)*(nImageN),(nIterN)*(nImageN-1));
		ptable.show("Results");
		ptableCC.show("Average CC");
		if(!bMultiCh)
		{
			MiscUtils.wrapFloatImgCal(averageImg,"final average "+Integer.toString(iter+1),cal,false).show();
		}
		else
		{
			showMultiChAverage(imgs_multiCh,shifts,"final average "+Integer.toString(iter+1),cal);
		}
	}
	
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
	
	Img< FloatType > subtractFromAverage(RandomAccessibleInterval< FloatType > aver_image, RandomAccessibleInterval< FloatType > curr_image)
	{
		int nDim = aver_image.numDimensions();
		// here is strong assumption that origin of coordinates is at 0
		final Img<FloatType> finImg = ArrayImgs.floats(aver_image.dimensionsAsLongArray());
		
		RandomAccessibleInterval< FloatType > subImg = Views.interval(Views.extendZero(curr_image), aver_image);
		RandomAccess<FloatType> raAver = aver_image.randomAccess();
		RandomAccess<FloatType> raSub = subImg.randomAccess();
		final Cursor< FloatType > finC = finImg.cursor();
		long [] pos = new long [nDim]; 
		while(finC.hasNext())
		{
			finC.fwd();
			finC.localize(pos);
			raAver.setPosition(pos);
			raSub.setPosition(pos);
			finC.get().set(raAver.get().get()-raSub.get().get());
		}
		
		return finImg;
	}
}
