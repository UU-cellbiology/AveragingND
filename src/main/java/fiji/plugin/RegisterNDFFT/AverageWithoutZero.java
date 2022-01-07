package fiji.plugin.RegisterNDFFT;

import java.util.ArrayList;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.plugin.PlugIn;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.ImagePlusAdapter;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;

public class AverageWithoutZero implements PlugIn {

	public boolean multiCh=false;
	@Override
	public void run(String arg0) {

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

		ArrayList<RandomAccessibleInterval< FloatType >> imgs = new ArrayList<RandomAccessibleInterval< FloatType >>();
		
		for(i=0;i<idList.length;i++)
		{
			imgs.add(ImagePlusAdapter.convertFloat(WindowManager.getImage(idList[i]))); 
		}

		//ImagePlus averagedIP = ImageJFunctions.wrap(averageImg,"Average");
		IntervalView<FloatType> averageImg = averageArray(imgs);
		if(WindowManager.getImage(0).getNChannels()>1)
		{
			multiCh = true;
		}
		MiscUtils.wrapFloatImgCal(averageImg,"Averaged_no_zeros", WindowManager.getImage(idList[0]).getCalibration(),multiCh).show();
		/*ImagePlus averagedIP = ImageJFunctions.wrap(averageImg,"Averaged_no_zeros");
		//ImagePlus registeredT = ImageJFunctions.wrap(Views.interval(Views.translate(Views.extendValue(template,Float.NaN), shift_pos),image),"GNCC registered_template");
		if(averageImg.numDimensions()==3)
		{
			averagedIP.setDimensions(1, (int)averageImg.dimension(2), 1);
		}
		if(averageImg.numDimensions()==4)
		{
			averagedIP.setDimensions((int)averageImg.dimension(2), (int)averageImg.dimension(3), 1);
		}
		averagedIP.show();
		*/
	}
	public static FinalInterval getIntervalAverageArray(ArrayList<RandomAccessibleInterval< FloatType >> imgs)
	{
		int i,j;
		
		int nDim = imgs.get(0).numDimensions();
		long [][] maxDim = new long[2][nDim];
		long [] currMinDim = new long[nDim];
		long [] currMaxDim = new long[nDim];
		for(i=0;i<imgs.size();i++)
		{
			imgs.get(i).max(currMaxDim);
			imgs.get(i).min(currMinDim);
			//imgs.get(i).dimensions(currDim);
			for (j=0;j<nDim;j++)
			{
				if(maxDim[0][j]>currMinDim[j])
				{
					maxDim[0][j]=currMinDim[j];
				}

				if(maxDim[1][j]<currMaxDim[j])
				{
					maxDim[1][j]=currMaxDim[j];
				}
			}
		}
		long [] dimMax = new long[nDim];
		for(j=0;j<nDim;j++)
		{
			dimMax[j]=maxDim[1][j]-maxDim[0][j]+1;
		}
		return new FinalInterval( maxDim[0] ,  maxDim[1] );
	}
	
	
	public static IntervalView<FloatType> averageArray(ArrayList<RandomAccessibleInterval< FloatType >> imgs)
	{
		int i;
		int nDim = imgs.get(0).numDimensions();

		FinalInterval intervalMax = getIntervalAverageArray(imgs);
		
		ArrayList<IntervalView< FloatType >> interv = new ArrayList<IntervalView< FloatType >>();
		for(i=0;i<imgs.size();i++)
		{
			interv.add(Views.interval( Views.extendZero(imgs.get(i)),intervalMax));
		}
		
		ArrayList<Cursor< FloatType >> cursors = new ArrayList<Cursor< FloatType >>();
		for(i=0;i<interv.size();i++)
		{
			cursors.add(interv.get(i).cursor());
		}
		final Img<FloatType> averageImgArr = ArrayImgs.floats(intervalMax.dimensionsAsLongArray());
		long [] originCoord = intervalMax.minAsLongArray();
		final IntervalView<FloatType> averageImg = Views.translate(averageImgArr, originCoord);
		//final ArrayImg<FloatType, FloatArray> averageImg = ArrayImgs.floats(dimMax);
		Cursor<FloatType> avC = averageImg.cursor();
		Cursor<FloatType> imgC;
		float nNum, nVal, nValCur;
		while(avC.hasNext())
		{
			avC.fwd();
			nNum=0.0f;
			nVal=0.0f;
			for(i=0;i<cursors.size();i++)
			{
				imgC=cursors.get(i);
				imgC.fwd();
				nValCur=imgC.get().get();
				if(nValCur>0)
				{
					nNum++;
					nVal+=nValCur;
				}
			}
			if(nNum>0.0f)
			{
				avC.get().set(nVal/nNum);
			}
		}	
		return 	Views.zeroMin(averageImg);

	}
	

	
	public static ArrayList<IntervalView<FloatType>> sumAndCountArray(ArrayList<RandomAccessibleInterval< FloatType >> imgs)
	{
		int i;
		//int nDim = imgs.get(0).numDimensions();

		FinalInterval intervalMax = getIntervalAverageArray(imgs);
		
		ArrayList<IntervalView< FloatType >> interv = new ArrayList<IntervalView< FloatType >>();
		for(i=0;i<imgs.size();i++)
		{
			interv.add(Views.interval( Views.extendZero(imgs.get(i)),intervalMax));
		}
		
		ArrayList<Cursor< FloatType >> cursors = new ArrayList<Cursor< FloatType >>();
		for(i=0;i<interv.size();i++)
		{
			cursors.add(interv.get(i).cursor());
		}
		final Img<FloatType> sumImgArr = ArrayImgs.floats(intervalMax.dimensionsAsLongArray());
		final Img<FloatType> countImgArr = ArrayImgs.floats(intervalMax.dimensionsAsLongArray());
		long [] originCoord = intervalMax.minAsLongArray();
		final IntervalView<FloatType> sumImg = Views.translate(sumImgArr, originCoord);
		final IntervalView<FloatType> countImg = Views.translate(countImgArr, originCoord);
		//final ArrayImg<FloatType, FloatArray> averageImg = ArrayImgs.floats(dimMax);
		Cursor<FloatType> sumC = sumImg.cursor();
		Cursor<FloatType> cntC = countImg.cursor();
		Cursor<FloatType> imgC;
		float nNum, nVal, nValCur;
		while(sumC.hasNext())
		{
			sumC.fwd();
			cntC.fwd();
			nNum=0.0f;
			nVal=0.0f;
			for(i=0;i<cursors.size();i++)
			{
				imgC=cursors.get(i);
				imgC.fwd();
				nValCur=imgC.get().get();
				if(nValCur>0)
				{
					nNum++;
					nVal+=nValCur;
				}
			}
			if(nNum>0.0f)
			{
				sumC.get().set(nVal);
				cntC.get().set(nNum);
			}
		}

		ArrayList<IntervalView<FloatType>> finalSumCnt = new ArrayList<IntervalView<FloatType>>();
		finalSumCnt.add(sumImg);
		finalSumCnt.add(countImg);
		return 	finalSumCnt;

	}
	
	public static IntervalView<FloatType> averageFromSumAndCount(ArrayList<IntervalView<FloatType>> alSumCnt)
	{

		long [] origin = alSumCnt.get(0).minAsLongArray();
		final Img<FloatType> avrgImgArr = ArrayImgs.floats(alSumCnt.get(0).dimensionsAsLongArray());
		final IntervalView<FloatType> avrgImg = Views.translate(avrgImgArr, origin );
		Cursor<FloatType> avrgC = avrgImg.cursor();
		Cursor<FloatType> sumC = alSumCnt.get(0).cursor();
		Cursor<FloatType> cntC = alSumCnt.get(1).cursor();
		float fCnt;
		while(avrgC.hasNext())
		{
			avrgC.fwd();
			sumC.fwd();
			cntC.fwd();
			fCnt=cntC.get().get();
			if(fCnt>0.0f)
			{
				avrgC.get().set(sumC.get().get()/fCnt);
			}
		}				
		
		return Views.zeroMin(avrgImg);
	}
}
