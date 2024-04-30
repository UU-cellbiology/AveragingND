package averagingND;



import java.awt.AWTEvent;
import java.awt.Choice;
import java.awt.Label;
import java.awt.TextField;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;

import ij.IJ;
import ij.ImageJ;
import ij.Prefs;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.io.DirectoryChooser;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import net.imglib2.FinalInterval;


public class PairWiseCC implements PlugIn, DialogListener {

	public int nInput = 0; 

	boolean bZeroMask = false;
		
	/** set of images for averaging and information about them **/
	ImageSet imageSet;
	
	public int nConstrainReg = 0;
	
	/** choice UI for constrain type **/
	Choice limitCh;

	/** labels of constrain axes **/
	Label [] limName;
	
	/** values of constrain axes **/
	TextField [] limVal;
	
	int nDimReg;
	String sDims;
	boolean bCenteredLimit = true;
	
	@Override
	public void run(String arg) {
		
		
		int i,j,k;
		int d;
		
		//double format formatting tool
		DecimalFormatSymbols symbols = new DecimalFormatSymbols();
		symbols.setDecimalSeparator('.');
		//DecimalFormat df = new DecimalFormat ("#.########", symbols);
		DecimalFormat df1 = new DecimalFormat ("#.#", symbols);
		
		double [] dLimits;
		final String[] sInput = new String[2];
		sInput[0] = "All currently open images";
		sInput[1] = "Specify images in a folder";
		
		final GenericDialog gdFiles = new GenericDialog( "Pairwise CC" );
		gdFiles.addChoice( "Input images:", sInput, Prefs.get("RegisterNDFFT.PW.nInput", sInput[0]) );

		gdFiles.showDialog();
		
		if ( gdFiles.wasCanceled() )
			return;		
		IJ.log("Iterative ND averaging plugin, version " + ConstantsAveragingND.sVersion);
		IJ.log("Pairwise CC (cross-correlation) command.");
		
		nInput = gdFiles.getNextChoiceIndex();
		Prefs.set("RegisterNDFFT.PW.nInput", sInput[nInput]);


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
		
		sDims = imageSet.sRefDims;
		
		nDimReg = sDims.length();
		if(imageSet.bMultiCh)
		{
			nDimReg--; //remove the C component
		}
		limName = new Label[nDimReg];
		limVal = new TextField[nDimReg];
		dLimits = new double [nDimReg];
		
		final String[] limitsReg = new String[  ] {"No","by voxels", "by image fraction"};
		final GenericDialog gd1 = new GenericDialog( "CC parameters" );
		if(imageSet.bMultiCh)
		{
			final String[] channels = new String[ imageSet.nChannels];
			for ( int c = 0; c < channels.length; ++c )
				channels[ c ] = "use channel " + Integer.toString(c+1);
			gd1.addChoice( "For calculations ", channels, channels[ 0 ] );
		}
		
		gd1.addCheckbox("Use zero masked CC?", Prefs.get("RegisterNDFFT.PW.bExcludeZeros", false));	
		String sCurrChoice = Prefs.get("RegisterNDFFT.PW.sConstrain", "No");
		gd1.addChoice("Constrain registration?", limitsReg, sCurrChoice);
		limitCh = (Choice) gd1.getChoices().lastElement();
		for (d=0;d<nDimReg;d++)
		{
			switch (sCurrChoice)
			{
				case "No":
					gd1.addNumericField("No max "+sDims.charAt(d)+" limit", 0.0, 3);
					break;
				case "by voxels":
					gd1.addNumericField(sDims.charAt(d)+" limit (px)", Prefs.get("RegisterNDFFT.PW.dMax"+sDims.charAt(d)+"px", 10.0), 3);
					break;
				case "by image fraction":
					gd1.addNumericField(sDims.charAt(d)+" limit (0-1)", Prefs.get("RegisterNDFFT.PW.dMax"+sDims.charAt(d)+"fr", 0.5), 3);
					break;
					
			}
			limName[d] = gd1.getLabel();
			limVal[d] = (TextField)gd1.getNumericFields().get(d);	
			if(sCurrChoice.equals("No"))
			{
				limVal[d].setEnabled(false);
			}
		}
		
		gd1.addCheckbox("Image centered constrains?", Prefs.get("RegisterNDFFT.PW.bCenteredLimit", false));

		gd1.addDialogListener(this);
		gd1.showDialog();
		
		if ( gd1.wasCanceled() )
			return;

		
		if(imageSet.bMultiCh)
		{
			imageSet.alignChannel = gd1.getNextChoiceIndex();
		}

		bZeroMask  = gd1.getNextBoolean();
		Prefs.set("RegisterNDFFT.PW.bExcludeZeros", bZeroMask);
		
		nConstrainReg = gd1.getNextChoiceIndex();
		Prefs.set("RegisterNDFFT.PW.sConstrain", limitsReg[nConstrainReg]);

		
		bCenteredLimit  = gd1.getNextBoolean();
		Prefs.set("RegisterNDFFT.PW.bCenteredLimit", bCenteredLimit);

		if(nConstrainReg==0)
		{
			IJ.log("Pairwise CC without constrains.");
		}
		else		
		{
			if(nConstrainReg == 1)
			{
				IJ.log("Pairwice CC with constrain specified in voxels:");
				
				for(d=0;d<nDimReg;d++)
				{
					dLimits[d]=Math.abs(gd1.getNextNumber());
					Prefs.set("RegisterNDFFT.PW.dMax"+sDims.charAt(d)+"px",dLimits[d]);
					IJ.log("Axis " +sDims.charAt(d)+": "+df1.format(dLimits[d])+" pixels");
				}
				
			}
			else
			{
				IJ.log("Averaging with constrain specified as a fraction of max displacement:");
				for(d=0;d<nDimReg;d++)
				{
					dLimits[d]=Math.min(Math.abs(gd1.getNextNumber()), 1.0);
					Prefs.set("RegisterNDFFT.PW.dMax"+sDims.charAt(d)+"fr",dLimits[d]);
					IJ.log("Axis " +sDims.charAt(d)+": "+df1.format(dLimits[d]));
				}
			}
			if(bCenteredLimit)
			{
				IJ.log("Constrains applied with respect to centered position.");
			}
			else
			{
				IJ.log("Constrains applied with respect to the coordinates origin (Zero, top-left).");
			}
		}
		
		double [] lim_fractions = null;
		FinalInterval limInterval = null;
		if(nConstrainReg == 1)
		{
			long[] minI = new long [nDimReg];
			long[] maxI = new long [nDimReg];
			for(d=0;d<nDimReg;d++)
			{
				maxI[d] = (long) dLimits[d];
				minI[d] = (long) ((-1.0)*dLimits[d]);
			}
			limInterval = new FinalInterval(minI, maxI);
		}
		if(nConstrainReg == 2)
		{
			lim_fractions = new double [nDimReg];
			for(d=0;d<nDimReg;d++)
			{
				lim_fractions[d] = dLimits[d];
			}
		}
		
		if(!imageSet.loadAllImages())
			return;
		
		final int nImageN = imageSet.nImageN;
		
		GenNormCC normCC = new GenNormCC();
		normCC.bVerbose = false;
		normCC.bZeroMask = bZeroMask;
		normCC.lim_fractions = lim_fractions;
		normCC.limInterval = limInterval;
		normCC.bCenteredLimit = bCenteredLimit;
		
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
				normCC.caclulateGenNormCC(imageSet.imgs.get(i), imageSet.imgs.get(j), false);
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

	@Override
	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
		int d;
		
		if(e!=null)
		{
			DecimalFormatSymbols symbols = new DecimalFormatSymbols();
			symbols.setDecimalSeparator('.');
			DecimalFormat df1 = new DecimalFormat ("#.##", symbols);

			if(e.getSource()==limitCh)
			{
				switch (limitCh.getSelectedIndex())
				{
					case 0:
						for(d=0;d<nDimReg;d++)
						{
							limName[d].setText("No "+sDims.charAt(d)+" limit");
							limVal[d].setEnabled(false);
						}
						break;
					case 1:
						for(d=0;d<nDimReg;d++)
						{
							limName[d].setText(sDims.charAt(d)+" limit (px)");
							limVal[d].setEnabled(true);
							limVal[d].setText(df1.format(Prefs.get("RegisterNDFFT.PW.dMax"+sDims.charAt(d)+"px", 10.0)));
						}
						break;
					case 2:
						for(d=0;d<nDimReg;d++)
						{
							limName[d].setText(sDims.charAt(d)+" limit (0-1)");
							limVal[d].setEnabled(true);
							limVal[d].setText(df1.format(Prefs.get("RegisterNDFFT.PW.dMax"+sDims.charAt(d)+"fr", 0.5)));

						}
						break;
				}
			}
		}
		return true;
	}
	
	public static void main( final String[] args )
	{
		// open an ImageJ window
		 new ImageJ();
		//IJ.open("/home/eugene/Desktop/projects/RegisterNDFFT/single/MAX_089-1.tif");
		//IJ.open("/home/eugene/Desktop/projects/RegisterNDFFT/single/MAX_098-1.tif");
		 PairWiseCC pw = new PairWiseCC();
		 pw.run(null);

	}

	

}
