package averagingND;

import ij.ImagePlus;
import net.imglib2.cache.img.CellLoader;
import net.imglib2.cache.img.ReadOnlyCachedCellImgFactory;
import net.imglib2.cache.img.ReadOnlyCachedCellImgOptions;
import net.imglib2.cache.img.SingleCellArrayImg;
import net.imglib2.img.Img;
import net.imglib2.type.numeric.real.FloatType;

public class FloatTiffImgWrap {

	
	public static Img< FloatType > wrapVirtualFloat(final ImagePlus imp, String sDims_)
	{
		String sDims = sDims_;
		final int nDims = sDims.length();
		final long[] dimensions = new long[nDims];
		final int nBitD = imp.getBitDepth();
		for(int i=0;i<nDims;i++)
		{
			if(sDims.charAt(i)=='X')
			{
				dimensions[i]=imp.getWidth();
			}
			if(sDims.charAt(i)=='Y')
			{
				dimensions[i]=imp.getHeight();
			}
			if(sDims.charAt(i)=='Z')
			{
				dimensions[i]=imp.getNSlices();
			}
			if(sDims.charAt(i)=='C')
			{
				dimensions[i]=imp.getNChannels();				
			}
			if(sDims.charAt(i)=='T')
			{
				dimensions[i]=imp.getNFrames();
			}
		}
		
		// set up cell size such that one cell is one plane
		final int[] cellDimensions = new int[] {
				imp.getStack().getWidth(),
				imp.getStack().getHeight(),
				1
		};
		// make a CellLoader that copies one plane of data from the virtual stack
		final CellLoader< FloatType > loader = new CellLoader< FloatType  >()
		{
			@Override
			public void load( final SingleCellArrayImg< FloatType , ? > cell ) throws Exception
			{
				//final int z = ( int ) cell.min( 2 );
				int nCh = 1;
				int nSl = 1;
				int nTp = 1;
				for(int i=2;i<nDims;i++)
				{
					if(sDims.charAt(i)=='C')
					{
						nCh = (int)cell.min(i)+1;
					}
					if(sDims.charAt(i)=='Z')
					{
						nSl = (int)cell.min(i)+1;
					}
					if(sDims.charAt(i)=='T')
					{
						nTp = (int)cell.min(i)+1;
					}

				}
				
				final int stInd = imp.getStackIndex(nCh, nSl, nTp);
				//final short[] impdata = ( short[] ) imp.getStack().getProcessor( 1 + z ).getPixels();
				final float[] celldata = ( float[] ) cell.getStorageArray();
				if(nBitD == 32)
				{
					final float[] impdata = ( float[] ) imp.getStack().getProcessor(stInd).getPixels();
					System.arraycopy( impdata, 0, celldata, 0, celldata.length );
				}
				if(nBitD == 16)
				{
					final short[] impdata = ( short[] ) imp.getStack().getProcessor( stInd ).getPixels();
					for(int i=0;i<celldata.length;i++)
					{
						celldata[i]=impdata[i]&0xffff;
					}
				}
				if(nBitD == 8)
				{
					final byte[] impdata = ( byte[] ) imp.getStack().getProcessor( stInd ).getPixels();
					for(int i=0;i<celldata.length;i++)
					{
						celldata[i]=impdata[i]&0xff;
					}
				}
				//System.arraycopy( impdata, 0, celldata, 0, celldata.length );
			}
		};
		// create a CellImg with that CellLoader
//		final Img< FloatType > img = new ReadOnlyCachedCellImgFactory().create(
//				dimensions,
//				new FloatType(),
//				loader,
//				ReadOnlyCachedCellImgOptions.options().cellDimensions( cellDimensions ) );
		return new ReadOnlyCachedCellImgFactory().create(
				dimensions,
				new FloatType(),
				loader,
				ReadOnlyCachedCellImgOptions.options().cellDimensions( cellDimensions ) );
	}
	//static long[]
}
