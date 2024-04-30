package averagingND;

import java.util.ArrayList;

import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Intervals;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;

public class TemplateAveraging {
	
	/** array holding information about current template.
	 * In case of average, contains a sum of all pixels,
	 * in case of median - three imgs around the median
	 * in case of masked average - sum and count per pixel **/
	public ArrayList<IntervalView<FloatType>> current_stat = new ArrayList<IntervalView<FloatType>>();
	
	
	public IntervalView<FloatType> currentTemplate = null;
	
	/** possible types of ImageSet statistics stored **/
	public static final int AVERAGE=0, MASKED_AVERAGE=1, MEDIAN=2; 
	
	/** specific type that is used to analyze data**/
	public int nStatType;
	
	/** if the object was initialized **/
	public boolean bInit = false;
	
	public int nImgN = 0;
	
	public FinalInterval unionInterval;
	
	public TemplateAveraging(final int nStatType_)
	{
		nStatType = nStatType_;
	}
	
	/** initialize the object **/
	public void init(ArrayList<RandomAccessibleInterval< FloatType >> imgs_)
	{
		switch(nStatType)
		{
			case AVERAGE:
				unionInterval = sumArray(imgs_, current_stat);
				break;
			case MASKED_AVERAGE:
				unionInterval = sumAndCountArray(imgs_, current_stat);
				break;
			case MEDIAN:
				break;
		}
		currentTemplate = null;
		bInit = true;
		nImgN = imgs_.size();
	}
	
	public IntervalView<FloatType> getTemplateForImage(RandomAccessibleInterval< FloatType > currentImage)
	{
		IntervalView<FloatType> out = null;

		if(!bInit)
		{
			System.out.println("Template object was not initialized!");
			return out;
		}
		
		switch(nStatType)
		{
			case AVERAGE:
				out = getTemplateForImageAverage(currentImage);
				break;
			case MASKED_AVERAGE:
				out = getTemplateForImageMaskedAverage(currentImage);
				break;
			case MEDIAN:
				break;
		}
		
		return out;
	}
	
	/** allocates memory for the current template **/
	void initCurrentTemplate()
	{
		final long [] origin = current_stat.get(0).minAsLongArray();
		final Img<FloatType> avrgImgArr = ArrayImgs.floats(current_stat.get(0).dimensionsAsLongArray());
		currentTemplate = Views.translate(avrgImgArr, origin );
	}
	
	/** returns current template for AVERAGE,
	 *  i.e. sum - currentImage/number of images **/
	IntervalView<FloatType> getTemplateForImageAverage(final RandomAccessibleInterval< FloatType > currentImage)
	{
	
		if(currentTemplate == null)
		{
			initCurrentTemplate();
		}
		final float nIm = nImgN-1.0f;
		
		final IntervalView<FloatType> removeInt = Views.interval(Views.extendZero(currentImage), currentTemplate);
		final Cursor<FloatType> avrgC = currentTemplate.cursor();
		final Cursor<FloatType> remC = removeInt.cursor();
		final Cursor<FloatType> sumC = current_stat.get(0).cursor();

		while(avrgC.hasNext())
		{			
			avrgC.fwd();
			remC.fwd();
			sumC.fwd();
			avrgC.get().set((sumC.get().get()-remC.get().get())/nIm);
		}
		//ImageJFunctions.show(currentTemplate, "Xffs");	
		//return Views.zeroMin(avrgImg);
		return currentTemplate;
	}
	
	/** returns current template for MASKED_AVERAGE,
	 *  i.e. sum - currentImage/updated count per pixel **/
	IntervalView<FloatType> getTemplateForImageMaskedAverage(final RandomAccessibleInterval< FloatType > currentImage)
	{
	
		if(currentTemplate == null)
		{
			initCurrentTemplate();
		}
		
		final IntervalView<FloatType> removeInt = Views.interval(Views.extendZero(currentImage), currentTemplate);
		final Cursor<FloatType> avrgC = currentTemplate.cursor();
		final Cursor<FloatType> remC = removeInt.cursor();
		final Cursor<FloatType> sumC = current_stat.get(0).cursor();
		final Cursor<FloatType> cntC = current_stat.get(1).cursor();
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
		return currentTemplate;
	}

	
	/** Provided with an ArrayList of RAIs, makes a new current_stat ArrayList with a single RAI,
	 * containing cumulative sum of all intensities at a current voxel/pixel location. 
	 * Since input RAIs could be of different sizes/locations, the output size is made to include
	 * all of them, i.e. a hyper box that includes them all (stored in unionInterval) **/	
	public static FinalInterval sumArray(final ArrayList<RandomAccessibleInterval< FloatType >> imgs, final ArrayList<IntervalView< FloatType >> out)
	{
		int i;

		final FinalInterval outUnionInterval = getUnionIntervalFromArray(imgs);
		
		final ArrayList<IntervalView< FloatType >> interv = new ArrayList<IntervalView< FloatType >>();
		
		for(i=0;i<imgs.size();i++)
		{
			interv.add(Views.interval( Views.extendZero(imgs.get(i)),outUnionInterval));
		}
		
		ArrayList<Cursor< FloatType >> cursors = new ArrayList<Cursor< FloatType >>();
		for(i=0;i<interv.size();i++)
		{
			cursors.add(interv.get(i).cursor());
		}
		
		final Img<FloatType> sumImgArr = ArrayImgs.floats(outUnionInterval.dimensionsAsLongArray());
		
		
		long [] originCoord = outUnionInterval.minAsLongArray();
		
		final IntervalView<FloatType> sumImg = Views.translate(sumImgArr, originCoord);


		final Cursor<FloatType> sumC = sumImg.cursor();	
		Cursor<FloatType> imgC;
		double nSumVal;		
		while(sumC.hasNext())
		{
			sumC.fwd();
			nSumVal = 0;
			for(i=0;i<cursors.size();i++)
			{
				imgC = cursors.get(i);
				imgC.fwd();
				nSumVal += imgC.get().get();
			}
			sumC.get().set((float)nSumVal);

		}
		out.clear();
		out.add(sumImg);
		return outUnionInterval;

	}
	
	/** Provided with an ArrayList of RAIs, makes a new current_stat ArrayList with two RAIs:
	 * the first contains cumulative sum of all intensities at a current voxel/pixel location,
	 * the second contains an integer value equal to how many RAIs have a pixel at this location. 
	 * Since input RAIs could be of different sizes/locations, the output size is made to include
	 * all of them, i.e. a hyper box that includes them all (stored in unionInterval) **/
	
	public static FinalInterval sumAndCountArray(ArrayList<RandomAccessibleInterval< FloatType >> imgs, final ArrayList<IntervalView< FloatType >> out)
	{
		int i;

		final FinalInterval outUnionInterval = getUnionIntervalFromArray(imgs);
		
		
		ArrayList<IntervalView< FloatType >> interv = new ArrayList<IntervalView< FloatType >>();
		
		for(i=0;i<imgs.size();i++)
		{
			interv.add(Views.interval( Views.extendZero(imgs.get(i)),outUnionInterval));
		}
		
		ArrayList<Cursor< FloatType >> cursors = new ArrayList<Cursor< FloatType >>();
		for(i=0;i<interv.size();i++)
		{
			cursors.add(interv.get(i).cursor());
		}
		
		final Img<FloatType> sumImgArr = ArrayImgs.floats(outUnionInterval.dimensionsAsLongArray());
		final Img<FloatType> countImgArr = ArrayImgs.floats(outUnionInterval.dimensionsAsLongArray());
		
		long [] originCoord = outUnionInterval.minAsLongArray();
		
		final IntervalView<FloatType> sumImg = Views.translate(sumImgArr, originCoord);
		final IntervalView<FloatType> countImg = Views.translate(countImgArr, originCoord);

		Cursor<FloatType> sumC = sumImg.cursor();
		Cursor<FloatType> cntC = countImg.cursor();
		Cursor<FloatType> imgC;
		long nNumVal;
		double nSumVal;
		float nValCur;
		
		while(sumC.hasNext())
		{
			sumC.fwd();
			cntC.fwd();
			nNumVal = 0;
			nSumVal = 0;
			for(i=0;i<cursors.size();i++)
			{
				imgC=cursors.get(i);
				imgC.fwd();
				nValCur=imgC.get().get();
				if(nValCur > 0.0000001)
				{
					nNumVal++;
					nSumVal += nValCur;
				}
			}
			if(nNumVal > 0)
			{
				sumC.get().set((float)nSumVal);
				cntC.get().set((float)nNumVal);
			}
		}

		out.clear();
		out.add(sumImg);
		out.add(countImg);
		return outUnionInterval;

	}
	
	/** returns an interval that encompasses all RAIs in the input ArrayList **/
	public static FinalInterval getUnionIntervalFromArray(ArrayList<RandomAccessibleInterval< FloatType >> imgs)
	{
		
		if (imgs.size()==0)
			return null;
		FinalInterval out = new FinalInterval(imgs.get(0));
		
		for(int i=1;i<imgs.size();i++)
			out = Intervals.union(out, imgs.get(i));

		return out;
	}
}
