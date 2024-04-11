package registerNDFFT;



import ij.IJ;

import ij.Prefs;

import ij.gui.GenericDialog;
import ij.io.DirectoryChooser;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;


public class PairWiseCC implements PlugIn {

	public int nInput = 0; 

	boolean bExcludeZeros = false;

	public double dMaxFraction = 0.5;
		
	/** set of images for averaging and information about them **/
	ImageSet imageSet;
	
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
		imageSet = new ImageSet();

		if(nInput == 0)
		{
			if(!imageSet.initializeFromOpenWindows())
				return;	 
		}
		else
		{
			DirectoryChooser dc = new DirectoryChooser ( "Choose a folder with images.." );
			String sPath = dc.getDirectory();
			if(!imageSet.initializeFromDisk(sPath, ".tif"))
				return;
		}
		if(imageSet.bMultiCh)
		{
			final GenericDialog gdCh = new GenericDialog( "Choose registration channel" );
		
			final String[] channels = new String[ imageSet.nChannels];
			for ( int c = 0; c < channels.length; ++c )
				channels[ c ] = "use channel " + Integer.toString(c+1);
			gdCh.addChoice( "For alignment", channels, channels[ 0 ] );
			gdCh.showDialog();
			
			if ( gdCh.wasCanceled() )
				return;				
	
			imageSet.alignChannel = gdCh.getNextChoiceIndex();
		}
		
		if(!imageSet.loadAllImages())
			return;
		
		final int nImageN = imageSet.nImageN;
		
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
				normCC.caclulateGenNormCC(imageSet.imgs.get(i),imageSet.imgs.get(j), dMaxFraction , false);
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
			ptableFN.addValue("filename", imageSet.image_names.get(i));
		}
		ptableFN.incrementCounter();
		ptableFN.addValue("filename", imageSet.image_names.get(i));
		
		IJ.showStatus("Calculating pairwise CC...done.");
		IJ.showProgress(2, 2);
		
		ptable.show("Results");
		ptableFN.show("Filenames");
		
	}
	

}
