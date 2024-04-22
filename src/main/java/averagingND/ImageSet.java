package averagingND;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.measure.Calibration;

import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.ImagePlusAdapter;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

/** class containing a collection of images with info about them,
 * plus function to load them **/
public class ImageSet {
	
	/** original images (with full channels) **/
	public ArrayList<RandomAccessibleInterval< FloatType >> imgs_multiCh;
	
	/** only channel used for the alignment **/
	public ArrayList<RandomAccessibleInterval< FloatType >> imgs;

	/** filenames of images **/
	public ArrayList<String> image_names;
	
	/** one or multiple channels images **/
	public boolean bMultiCh = false;
	
	/** total number of images **/
	public int nImageN = 0;
	
	public int nChannels;
	public int nSlices;
	int nTimePoints;

	public int alignChannel = 1;
	
	/** calibration of the images **/
	public Calibration cal;
	
	List<String> filenames;
	
	String [] sFileNamesShort;
		
	String sFileExtension;
	
	int[] idList;
	public int nDim;
	
	/** source of images: 0 - all images are open in ImageJ
	 * 	1 - images form disk **/
	public int nSource;
	
	/** text representation of dimensions in the format of XYZCT **/
	public String sRefDims;
	
	/** text representation of dimensions in the format of XYCZT **/
	public String sRefDimsIJ;

	public ImageSet()
	{
		imgs_multiCh = new ArrayList<RandomAccessibleInterval< FloatType >>();
		imgs = new ArrayList<RandomAccessibleInterval< FloatType >>();		
		image_names = new ArrayList<String>();
	}
	/** function reads first file from the disk, analyzes its dimensions
	 * and assumes the rest of files is the same**/			//if(!loadFolderTiff())
	//	return;
	boolean initializeFromDisk(String sPath, String sFileExtension_)
	{
		nSource = 1;
		sFileExtension = sFileExtension_;
		// getting list of files
		filenames = null;
		if (sPath != null)
		{
			IJ.log("Analyzing folder "+sPath);
			try {
				filenames = getFilenamesFromFolder(sPath, sFileExtension);
			} catch (IOException e) {
				e.printStackTrace();
				IJ.log(e.getMessage());
			}	
			if(filenames == null)
			{
				return false;
			}
			else
			if(filenames.isEmpty())
			{
				IJ.log("Cannot find any "+sFileExtension+" files in provided folder, aborting.");
				return false;
			}
			else
			{
				IJ.log("Found " +Integer.toString(filenames.size())+" "+sFileExtension+" files.");
			}
		}
		else
			return false;
		
		fillFilenamesArray(filenames);
		ImagePlus ipFirst = IJ.openImage(filenames.get(0));
		fillDimensions(ipFirst);
		ipFirst.close();

		return true;
	}
	
	
	boolean initializeFromOpenWindows()
	{
		nSource = 0;

		idList = WindowManager.getIDList();		

		if ( idList == null || idList.length < 2 )
		{
			IJ.error( "The plugin requires at least two open images." );
			return false;
		}
				
		fillDimensions( WindowManager.getImage(idList[0]));

		
		return true;
	}
	
	
	boolean loadAllImages()
	{
		int i;
		
		nImageN = 0;
		ImagePlus imageIn;
		int nTryImages;
	
		if(nSource == 0)
			nTryImages = idList.length;
		else
			nTryImages = filenames.size();

		for(i=0;i<nTryImages;i++)
		{
			if(nSource==0)
			{
				imageIn = WindowManager.getImage(idList[i]);
			}
			else
			{
				imageIn = IJ.openImage(filenames.get(i));
			}
			if(checkConsistency(imageIn))
			{
				if(!bMultiCh)
				{
					imgs.add(ImagePlusAdapter.convertFloat(imageIn));				
				}
				else
				{
					imgs_multiCh.add(ImagePlusAdapter.convertFloat(imageIn));
					imgs.add(Views.hyperSlice(imgs_multiCh.get(i),2,alignChannel));
				}

				image_names.add(imageIn.getTitle());
				nImageN++;
			}
			else
			{
				return false;
			}
		}	

		//calculate min max 
		long [] dimsMin = new long [nDim];
		long [] dimsMax = new long [nDim];
		long [] currDim = new long [nDim];
		if(!bMultiCh)
		{
			imgs.get(0).dimensions(dimsMin);
			imgs.get(0).dimensions(dimsMax);
		}
		else
		{
			imgs_multiCh.get(0).dimensions(dimsMin);
			imgs_multiCh.get(0).dimensions(dimsMax);
		}
		for(i=1;i<nImageN;i++)
		{
			if(!bMultiCh)
			{
				imgs.get(i).dimensions(currDim);
			}
			else
			{
				imgs_multiCh.get(i).dimensions(currDim);
			}
			for (int d=0;d<nDim;d++)
			{
				if(dimsMin[d]>currDim[d])
				{
					dimsMin[d] = currDim[d];
				}
				if(dimsMax[d]<currDim[d])
				{
					dimsMax[d] = currDim[d];
				}
			}
		}
		for(int d=0;d<nDim;d++)
		{
			IJ.log("Axis "+sRefDimsIJ.charAt(d)+" min: "+Long.toString(dimsMin[d])+" max: "+Long.toString(dimsMax[d]));
		}
		IJ.showProgress(2,2);
		IJ.showStatus("Loading images..done.");	
		if(nImageN<2)
		{
			IJ.error( "The plugin requires at least two images with the same dimensions. Aborting." );
			return false;
		}
		if(!bMultiCh)
		{
			IJ.log("Averaging "+ Integer.toString(nImageN) + " images.");
		}
		else
		{
			IJ.log("Averaging "+Integer.toString(nImageN) + " images with " + Integer.toString(nChannels) + " channels.");
			IJ.log("Using channel "+Integer.toString(alignChannel+1) + " for alignment.");
		}
		return true;
	}
	
	boolean checkConsistency(final ImagePlus imageIn)
	{
		if(!Objects.equals(sRefDims, MiscUtils.getDimensionsText(imageIn)))
		{
			IJ.log("Error! Image \""+imageIn.getTitle() +"\" has different dimensions, skipping.");
			IJ.log("("+MiscUtils.getDimensionsText(imageIn) +" vs assumed " + sRefDims + ")");
			return false;
		}
		Calibration calIn = imageIn.getCalibration();
		if(!cal.equals(calIn))
		{
			IJ.log("Warning: image \""+imageIn.getTitle() +"\" has different voxel size, including, but maybe it needs to be checked.");
						
		}
		return true;
	}
	
	/** fills info about image from provided ImagePlus and closes it**/
	void fillDimensions(final ImagePlus ipFirst)
	{
		IJ.log("Analyzing dimensions:");
		sRefDims = MiscUtils.getDimensionsText(ipFirst);
		sRefDimsIJ = MiscUtils.getDimensionsTextImageJ(ipFirst);
		nDim = sRefDims.length();
		cal = ipFirst.getCalibration();
		String sDims = "XY";
		
		nChannels = ipFirst.getNChannels();
		nSlices = ipFirst.getNSlices();
		nTimePoints = ipFirst.getNFrames();
		if(nSlices>1)
		{
			sDims = sDims + "Z";
		}
		
		if(nTimePoints>1)
		{
			sDims = sDims + "T";
		}
		
		if(nChannels>1)
		{
			bMultiCh = true;
			sDims = sDims + "C";
		}


		sDims = sDims +" and " + Integer.toString(ipFirst.getBitDepth())+"-bit";
		
		if(bMultiCh)
		{
			sDims = sDims +" with "+ Integer.toString(nChannels)+" channels";
		}
		IJ.log(" - Inferring general dimensions/pixel sizes from \""+ipFirst.getTitle()+"\"");
		IJ.log(" - Assuming all files are "+sDims+" ");
	}
	

	/**given the path to folder and file extension, returns List of filenames strings **/
	public static List<String> getFilenamesFromFolder(final String sFolderPath, final String fileExtension)    
			throws IOException {
		final Path path = Paths.get(sFolderPath);

		if (!Files.isDirectory(path)) {
			throw new IllegalArgumentException("Path must be a directory!");
		}

		List<String> result = null;

		try (Stream<Path> stream = Files.list(path)) {
			result = stream
					.filter(p -> !Files.isDirectory(p))
					// this is a path, not string,
					// this only test if path end with a certain path
					//.filter(p -> p.endsWith(fileExtension))
					// convert path to string first
					.map(p -> p.toString())
					.filter(f -> f.endsWith(fileExtension))
					.collect(Collectors.toList());
		} catch (IOException e) {
			e.printStackTrace();
		}
		Collections.sort(result);
		return result;
	}
	
	public void fillFilenamesArray(List<String> fullPath)
	{
		sFileNamesShort = new String[fullPath.size()];
		Path p;
		for(int i=0;i<fullPath.size();i++)
		{
			 p = Paths.get(fullPath.get(i));
			 sFileNamesShort[i] = p.getFileName().toString();
			 sFileNamesShort[i] = sFileNamesShort[i].substring(0, sFileNamesShort[i].length()-sFileExtension.length());
		}
		
	}
}
