package averagingND;

import java.util.ArrayList;

import ij.IJ;

import ij.WindowManager;
import ij.plugin.PlugIn;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.ImagePlusAdapter;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Intervals;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;

public class FinalOutput implements PlugIn {

	public boolean multiCh = false;
	@Override
	public void run(String arg0) {

		int i;
		
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
		IntervalView<FloatType> averageImg = Views.zeroMin(averageArray(imgs, true));
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
	/** returns an interval that encompasses all RAIs in the input ArrayList **/
	public static FinalInterval getIntervalAverageArray(ArrayList<RandomAccessibleInterval< FloatType >> imgs)
	{
		
		if (imgs.size()==0)
			return null;
		FinalInterval out = new FinalInterval(imgs.get(0));
		
		for(int i=1;i<imgs.size();i++)
			out = Intervals.union(out, imgs.get(i));

		return out;
	}
	
	/** function calculates average image from the provided ArrayList of RAI intervals.
	 * The output image is extended to cover all intervals. Intervals are extended with zeros
	 * to the same extent.
	 * If bIgnoreZeros is true, it does not include in the averaging pixel values below 0.0000001
	 * **/
	public static IntervalView<FloatType> averageArray(final ArrayList<RandomAccessibleInterval< FloatType >> imgs, final boolean bIgnoreZero)
	{
		int i;

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

		Cursor<FloatType> avC = averageImg.cursor();
		Cursor<FloatType> imgC;
		long nNumVal;
		double dSumVal;
		float fValCur;
		while(avC.hasNext())
		{
			avC.fwd();
			nNumVal = 0;
			dSumVal = 0.0;
			for(i=0;i<cursors.size();i++)
			{
				imgC = cursors.get(i);
				imgC.fwd();
				fValCur = imgC.get().get();
				if(bIgnoreZero)
				{
					if(fValCur > 0.0000001)
					{
						nNumVal++;
						dSumVal += fValCur;
					}
				}
				else
				{
					nNumVal++;
					dSumVal += fValCur;
				}		
			}
			if(nNumVal > 0)
			{
				avC.get().set((float)(dSumVal/(double)nNumVal));
			}
		}	
		return averageImg;

	}
	
	/** function calculates standard deviation image from the provided ArrayList of RAI intervals
	 * and an average image RAI. 
	 * The output image is extended to cover all intervals. Intervals are extended with zeros
	 * to the same extent.
	 * If bIgnoreZeros is true, it does not include in the averaging pixel values below 0.0000001
	 * **/
	public static IntervalView<FloatType> stdArray(final ArrayList<RandomAccessibleInterval< FloatType >> imgs, final RandomAccessibleInterval< FloatType > avrgRAI, final boolean bIgnoreZero)
	{
		int i;

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
		final Img<FloatType> stdImgArr = ArrayImgs.floats(intervalMax.dimensionsAsLongArray());
		long [] originCoord = intervalMax.minAsLongArray();
		final IntervalView<FloatType> stdImg = Views.translate(stdImgArr, originCoord);

		Cursor<FloatType> stdC = stdImg.cursor();
		//just in case, use extendZero
		Cursor<FloatType> avrgC = Views.interval(Views.extendZero(avrgRAI),intervalMax).cursor();

		Cursor<FloatType> imgC;
		long nNumVal;
		double dSumVal;
		double dAvrgVal;
		float fValCur;
		while(stdC.hasNext())
		{
			stdC.fwd();
			avrgC.fwd();
			nNumVal = 0;
			dSumVal = 0.0;
			dAvrgVal = avrgC.get().get();
			for(i=0;i<cursors.size();i++)
			{
				imgC = cursors.get(i);
				imgC.fwd();
				fValCur = imgC.get().get();
				if(bIgnoreZero)
				{
					if(fValCur > 0.0000001)
					{
						nNumVal++;
						fValCur-=dAvrgVal;
						dSumVal += fValCur*fValCur;
					}
				}
				else
				{
					nNumVal++;
					fValCur-=  dAvrgVal;
					dSumVal += fValCur*fValCur;
				}		
			}
			if(nNumVal > 0)
			{
				stdC.get().set((float)Math.sqrt(dSumVal/(double)nNumVal));
			}
		}	
		return stdImg;

	}

}
