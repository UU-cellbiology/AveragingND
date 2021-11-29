package fiji.plugin.RegisterNDFFT;


import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;

//import ij.IJ;
//import ij.ImageJ;

//import net.imagej.ImageJ;

import io.scif.img.ImgIOException;
import io.scif.img.ImgOpener;
import net.imglib2.exception.IncompatibleTypeException;
import net.imglib2.img.ImagePlusAdapter;
import net.imglib2.img.Img;
import net.imglib2.img.imageplus.FloatImagePlus;
import net.imglib2.type.numeric.real.FloatType;


/**
 * Perform template matching by finding maximum of normalized cross-correlation
 * using convolution in the Fourier domain.
 *
 * @author Eugene Katrukha
 * 
 */
public class RegisterNDFFT implements PlugIn 
{
	
	public static int defaultImg1 = 0;
	public static int defaultImg2 = 1;
	public static boolean bShowCC = true;
	public static double max_fraction = 0.5;

	@Override
	public void run(String arg) {

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
		for ( int i = 0; i < idList.length; ++i )
			imgList[ i ] = WindowManager.getImage(idList[i]).getTitle();
		
		/**
		 * Dialog for choosing the images
		 */
		final GenericDialog gd1 = new GenericDialog( "Images/volumes registration" );
		
		if ( defaultImg1 >= imgList.length || defaultImg2 >= imgList.length )
		{
			defaultImg1 = 0;
			defaultImg2 = 1;
		}
		
		gd1.addChoice("First_image (reference)", imgList, imgList[ defaultImg1 ] );
		gd1.addChoice("Second_image (to register)", imgList, imgList[ defaultImg2 ] );
		gd1.addNumericField("Maximum shift (fraction, 0-1 range)", 0.9, 3);
		gd1.addCheckbox("Show cross-correlation", true);
				
		gd1.showDialog();
		
		if ( gd1.wasCanceled() )
			return;
		
		ImagePlus imp1 = WindowManager.getImage( idList[ defaultImg1 = gd1.getNextChoiceIndex() ] );		
		ImagePlus imp2 = WindowManager.getImage( idList[ defaultImg2 = gd1.getNextChoiceIndex() ] );	
		
		
		max_fraction = gd1.getNextNumber();
		bShowCC  = gd1.getNextBoolean();
		
		final Img< FloatType > image_in = ImagePlusAdapter.convertFloat(imp1);
		final Img< FloatType > template_in = ImagePlusAdapter.convertFloat(imp2);

		GenNormCC.caclulateGenNormCC(image_in, template_in, max_fraction, bShowCC);

		RegisterTranslation reg = new RegisterTranslation(image_in.numDimensions());
		
		reg.registerTranslation(image_in, template_in);
		
	}
	
	public static void main( final String[] args ) throws ImgIOException, IncompatibleTypeException
	{
		// open an ImageJ window
		new ImageJ();
		
		// open with SCIFIO ImgOpener as FloatTypes
		ImgOpener io = new ImgOpener();
		/*
		final Img< FloatType > image_in = io.openImgs( "Zstack 1.1-1-1.tif",
			new FloatType() ).get( 0 );
		//final Img< FloatType > template_in = io.openImgs( "Zstack 1.1-1-1.tif",
		final Img< FloatType > template_in = io.openImgs( "Zstack 1.2-2-1.tif",
			new FloatType() ).get( 0 );
		*/
			
		
		
		final Img< FloatType > image_in = io.openImgs( "linetest2Dn.tif",
				new FloatType() ).get( 0 );
		final Img< FloatType > template_in = io.openImgs( "linetest2Dn.tif",
				new FloatType() ).get( 0 );
		

/*	final Img< FloatType > image_in = io.openImgs( "linetest2Dn.tif",
				new FloatType() ).get( 0 );
		final Img< FloatType > template_in = io.openImgs( "linetest2Dn_shift.tif",
				new FloatType() ).get( 0 );
*/
		GenNormCC.caclulateGenNormCC(image_in, template_in, 1.0, true);
		//RegisterTranslation reg = new RegisterTranslation(image_in .numDimensions());
		// run the example
		//reg.registerTranslation(image_in, template_in);

		
	
	}


}
