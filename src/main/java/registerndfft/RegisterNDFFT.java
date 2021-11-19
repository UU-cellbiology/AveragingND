package registerndfft;

import ij.IJ;
import ij.ImageJ;

import io.scif.img.ImgIOException;
import io.scif.img.ImgOpener;


import net.imglib2.exception.IncompatibleTypeException;
import net.imglib2.img.Img;
import net.imglib2.type.numeric.real.FloatType;


/**
 * Perform template matching by finding maximum of normalized cross-correlation
 * using convolution in the Fourier domain.
 *
 * @author Eugene Katrukha
 * 
 */
public class RegisterNDFFT
{

	public static void main( final String[] args ) throws ImgIOException, IncompatibleTypeException
	{
		// open an ImageJ window
		new ImageJ();
		// open with SCIFIO ImgOpener as FloatTypes
		ImgOpener io = new ImgOpener();
		
		final Img< FloatType > image_in = io.openImgs( "Zstack 1.1-1-1.tif",
			new FloatType() ).get( 0 );
		final Img< FloatType > template_in = io.openImgs( "Zstack 1.2-2-1.tif",
			new FloatType() ).get( 0 );
		/*
		final Img< FloatType > image_in = io.openImgs( "linetest2Dn.tif",
				new FloatType() ).get( 0 );
		final Img< FloatType > template_in = io.openImgs( "linetest2Dn_shift.tif",
				new FloatType() ).get( 0 );
		 */			
		RegisterTranslation reg = new RegisterTranslation(image_in .numDimensions());
		// run the example
		reg.registerTranslation(image_in, template_in);

		
	
	}
}
