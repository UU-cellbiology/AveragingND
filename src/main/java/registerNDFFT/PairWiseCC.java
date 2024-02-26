package registerNDFFT;

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
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.ImagePlusAdapter;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

public class PairWiseCC implements PlugIn {

	public int nInput = 0; 
	boolean bExcludeZeros = false;
	public boolean bMultiCh = false;
	public double dMaxFraction = 0.5;
	public int nImageN = 0;
	public int numChannels = 1;
	public int alignChannel = 1;
	public ArrayList<String> filenames;
	
	/** original images (with full channels) **/
	ArrayList<RandomAccessibleInterval< FloatType >> imgs_multiCh;
	/** only channel used for CC calculation **/
	ArrayList<RandomAccessibleInterval< FloatType >> imgs;
	@Override
	public void run(String arg) {
		int i,j,k;
		
		final String[] sInput = new String[2];
		sInput[0] = "All currently open images";
		sInput[1] = "Specify images in a folder";
		
		final GenericDialog gd = new GenericDialog( "Pairwise CC" );
		gd.addChoice( "Input images:", sInput, Prefs.get("RegisterNDFFT.PW.nInput", sInput[0]) );
		gd.addCheckbox("Exclude zero values?", Prefs.get("RegisterNDFFT.PW.bExcludeZeros", false));	
		gd.addNumericField("Maximum shift (fraction, 0-1 range)", Prefs.get("RegisterNDFFT.PW.dMaxFraction", 0.4), 3);
		gd.showDialog();
		
		if ( gd.wasCanceled() )
			return;		
		
		nInput = gd.getNextChoiceIndex();
		Prefs.set("RegisterNDFFT.PW.nInput", sInput[nInput]);
		bExcludeZeros  = gd.getNextBoolean();
		Prefs.set("RegisterNDFFT.PW.bExcludeZeros", bExcludeZeros);
		dMaxFraction  = gd.getNextNumber();
		Prefs.set("RegisterNDFFT.PW.dMaxFraction", dMaxFraction);
		//init image arrays		
		imgs_multiCh = new ArrayList<RandomAccessibleInterval< FloatType >>();
		imgs = new ArrayList<RandomAccessibleInterval< FloatType >>();	
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
		
		
		GenNormCC normCC = new GenNormCC();
		normCC.bVerbose = false;
		normCC.bExcludeZeros=bExcludeZeros;
		ResultsTable ptable = ResultsTable.getResultsTable();
		ptable.reset();
		ResultsTable ptableFN = new ResultsTable();
		int nProgress=0;
		IJ.showStatus("Calculating pairwise CC...");
		IJ.showProgress(nProgress, (int)(nImageN*(nImageN-1)*0.5));
		for(i=0;i<(nImageN-1);i++)
		{
			for(j=i+1;j<nImageN;j++)
			{
				normCC.caclulateGenNormCC(imgs.get(i),imgs.get(j), dMaxFraction , false);
				ptable.incrementCounter();
				
				//ptable.addValue("particle_pair", Integer.toString(i+1)+"_"+Integer.toString(j+1));
				ptable.addValue("norm_CC_coeff", normCC.dMaxCC);
				for(k=0;k<normCC.dShift.length;k++)
				{
					ptable.addValue("coord_"+Integer.toString(k+1), normCC.dShift[k]);
				}
				ptable.addValue("ind1", i+1);
				ptable.addValue("ind2", j+1);
				nProgress++;
				IJ.showProgress(nProgress, (int)(nImageN*(nImageN-1)*0.5));
			}
			ptableFN.incrementCounter();
			ptableFN.addValue("filename", filenames.get(i));
		}
		ptableFN.incrementCounter();
		ptableFN.addValue("filename", filenames.get(i));
		
		IJ.showStatus("Calculating pairwise CC...done.");
		IJ.showProgress(2, 2);
		
		ptable.show("Results");
		ptableFN.show("Filenames");
		
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

		filenames = new ArrayList<String>();
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
			filenames.add(WindowManager.getImage(idList[i]).getTitle());
		}
		if(!bMultiCh)
		{
			IJ.log("Pairwise CC for "+ Integer.toString(nImageN) + " images.");
		}
		else
		{
			IJ.log("Pairwise CC for "+Integer.toString(nImageN) + " images with " + Integer.toString(numChannels) + " channels.");
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
		filenames = new ArrayList<String>();
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
				filenames.add(files.get(i).substring(sPath.length()));
			
				if(!bMultiCh)
				{
					imgs.add(ImagePlusAdapter.convertFloat(impBridge));
				}
				else
				{
					imgs_multiCh.add(ImagePlusAdapter.convertFloat(impBridge));
					imgs.add(Views.hyperSlice(imgs_multiCh.get(i),2,alignChannel));
				}
				IJ.showProgress(i, nImageN-1);
			}
			IJ.showProgress(2,2);
			IJ.showStatus("Loading images..done.");	
			if(!bMultiCh)
			{
				IJ.log("Pairwise CC for "+ Integer.toString(nImageN) + " images.");
			}
			else
			{
				IJ.log("Pairwise CC for "+Integer.toString(nImageN) + " images with " + Integer.toString(numChannels) + " channels.");
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
