import java.awt.*;
import java.io.*;
import java.util.Arrays;
import ij.plugin.*;
import ij.process.*;
import ij.gui.*;

public class Fiber_Orientation extends PlugInFrame implements ActionListener {
	// Start defining constants for UI 
	public static final String THRESHOLD        = "Threshold";
	public static final String FFT_AND_ROT      = "FFT and Rotate";
	public static final String ORIGINAL_FFT     = "Orig. FFT";
	public static final String ERODE            = "Erode";
	public static final String DILATE           = "Dilate";	
	public static final String AUTOMATIZATED    = "Automatic";
	public static final String ANALYZE          = "Analize";
	public static final String BINARY_THRESHOLD = "Binary Thres.";
	public static final String DRAW_ANGLE       = "Draw Angles";
	// END defining constants for UI
	// Variables for User interface panels objects
	private Panel panel;
	private static Frame instance;
	private TextField threshold_value_tf;
	private TextField degree_separation_value_tf;
	public  ImagePlus ourImage;
	public  ImagePlus fftImage;
	public  Plot plot;
	// END Variables for User interface panels objects
	private int previousID;//storage the ID of the previousimage 
	//Constructor Method
	public Fiber_Orientation() {
		super("Fiber Orientation");//init the object
		if (instance!=null) {
			this.instance.toFront();//Frame to front
			return;//already initzialized,method can end
		}
		this.instance = this;// Set instance as global variable
		addKeyListener(IJ.getInstance());//listen button press
		//init self Layout for UI
		this.setLayout(new FlowLayout());
		// Init the panel to display
		panel = new Panel();
		panel.setLayout(new GridLayout(0, 2, 5, 5));
		//Add UI objects to the instance 
		//display textFields
		addTextFieldThreshold("Threshold");
		addTextFieldDegreeSeparation("Degree Sep.");
		//init cehck box
		this.addCheckBox(BINARY_THRESHOLD);
		this.addCheckBox(DRAW_ANGLE);
		//init buttons 
		this.addButton(AUTOMATIZATED,false);
		this.addButton(FFT_AND_ROT,false);
		this.addButton(THRESHOLD,true);
		this.addButton(ORIGINAL_FFT,true);
		this.addButton(ERODE,true);
		this.addButton(DILATE,true);
		this.addButton(ANALYZE,true);
		//Display the panel with the UI objects inside
		this.add(panel);
		this.pack();
		GUI.center(this);
		this.setVisible(true);
	}
	//Own  method to create and add buttons
	void addButton(String label,boolean disabled) {
		Button b = new Button(label);
		b.setEnabled(!disabled);
		b.addActionListener(this);
		b.addKeyListener(IJ.getInstance());
		panel.add(b);
	}
	//Own  method to create and add check box
	void addCheckBox(String label){
		panel.add(new Checkbox(label));
	}
	//Own  method to create and add the threshold text-field
	void addTextFieldThreshold(String label) {
		//170 is the default threshold
		this.threshold_value_tf = new TextField("170");
		panel.add(new Label(label));
		panel.add(this.threshold_value_tf);
	}
	//Own  method to create and add the Degree Separation
	void addTextFieldDegreeSeparation(String label){
		//1 degree  of separation is the default value
		this.degree_separation_value_tf = new TextField("1");
		panel.add(new Label(label));
		panel.add(this.degree_separation_value_tf);
	}
	//MEthod to enable/disable buttons with a respective label 
	void enableButton(String label, boolean enable){
		for(int i =0; i< this.panel.getComponentCount(); i++){
			Component component = this.panel.getComponent(i);
			if ( component.getClass().equals(Button.class) ){
				Button b = (Button) component;
				if(b.getLabel().equals(label) ){
					b.setEnabled(enable);
					i = this.panel.getComponentCount();
				}
			}
		}
	}
	//Return the value of a checkbox of the respective label
	boolean readCheckBox(String label){
		for(int i =0; i< this.panel.getComponentCount(); i++){
			Component component = this.panel.getComponent(i);
			if( component.getClass().equals(Checkbox.class) ){
				Checkbox cb = (Checkbox) component;
				if(cb.getLabel().equals(label) )
					return cb.getState();
			}
		}
		return false;
	}

	//Listener when a buytton is pressed
	public void actionPerformed(ActionEvent e) {
		//Open selected image
		ImagePlus imp = WindowManager.getCurrentImage();
		//Detect when there is not image open
		if (imp==null) {
			IJ.beep();
			IJ.showStatus("No image!!!");
			previousID = 0;
			return;
		}
		if (!imp.lock())
			{previousID = 0; return;}//still not finish previous thread

		int id = imp.getID();//get current image ID
		if (id!=previousID)
			imp.getProcessor().snapshot();//obtain the snapshot

		previousID   = id;
		String label = e.getActionCommand();//get current button selected

		if (label==null)
			return;
		//create a thread to run the process
		new Runner(label, imp,this);
	}
	//To listen when the windows is closed to clean objects.
	public void processWindowEvent(WindowEvent e) {
		super.processWindowEvent(e);
		if (e.getID()==WindowEvent.WINDOW_CLOSING) {
			instance = null;	
		}
	}
	// Thread "Runner" inner class
	class Runner extends Thread { // inner class
		private String command;
		private ImagePlus imp;
		private Fiber_Orientation caller;
		public  int threshold_value   = 0;
		public  int degree_separation = 0;
		//Instance method
		Runner(String command, ImagePlus imp,Fiber_Orientation callerL) {
			super(command);
			// Init variables to make them available to other methods
			this.command = command;
			this.imp     = imp;
			this.caller  = callerL;
			this.threshold_value   = Integer.parseInt(callerL.threshold_value_tf.getText());
			this.degree_separation = Integer.parseInt(callerL.degree_separation_value_tf.getText());
			//set thread priority
			setPriority(Math.max(getPriority()-2, MIN_PRIORITY));
			//start the thread.. the "run"- method is called
			start();
		}
		public void run() {
			try {//Call the implemented method were the algortihm is implemented
				runCommand(command, imp);
			} catch(OutOfMemoryError e) {
				IJ.outOfMemory(command);
				if (imp!=null) imp.unlock();
			} catch(Exception e) {
				CharArrayWriter caw = new CharArrayWriter();
				PrintWriter pw = new PrintWriter(caw);
				e.printStackTrace(pw);
				IJ.log(caw.toString());
				IJ.showStatus("");
				if (imp!=null) imp.unlock();
			}
		}
	    //Implemented thread
		void runCommand(String command, ImagePlus imp) {
			//Init the ROi to evaluate, and the time to measure
			ImageProcessor ip = imp.getProcessor();
			long startTime    = System.currentTimeMillis();
			Roi roi           = imp.getRoi();
			if (command.equals(THRESHOLD))
				{roi = null; ip.resetRoi();}
			ImageProcessor mask =  roi!=null?roi.getMask():null;
			if (command.equals(FFT_AND_ROT)){
				//The actual FFT and rotation method is executed
				ImagePlus fft_image  = this.getFftandRotate(imp, ip);
				//Storage the  Original FFT 
				this.caller.ourImage = fft_image;
				this.storeOriginalFft();
				//enable the buttons
				this.caller.enableButton(THRESHOLD, true);
				this.caller.enableButton(ORIGINAL_FFT, true);
				this.caller.enableButton(ERODE, true);
				this.caller.enableButton(DILATE, true);
				this.caller.enableButton(ANALYZE, true);
			}else if (command.equals(ORIGINAL_FFT)){
				 this.restore_original_fft();
			}else if (command.equals(THRESHOLD)){
				ImagePlus image_thresholded  = this.threshold(this.caller.ourImage);
			}else if (command.equals(ERODE)){
				this.erode();
			}else if (command.equals(DILATE)){
				this.dilate();
			}else if (command.equals(ANALYZE)){
        		this.reportResult();
			}else if (command.equals(AUTOMATIZATED)){
				// GET FFT 
				ImagePlus fft_image  = this.getFftandRotate(imp, ip);
				this.caller.ourImage = fft_image;
				this.storeOriginalFft();
				// Then BINARIZE 
				this.threshold(this.caller.ourImage);
				for(int ii=0;ii<3;ii++){//erode dilation process
					this.erode();
					this.dilation();
				}
				this.reportResult();//Start Analizing to obtain the OFD, and peak density
			}
			if (mask!=null) ip.reset(mask);
			imp.updateAndDraw();
			imp.unlock();
		}
		// Analyse the FFT spectrum
		void reportResult(){
			double[] resultsOFD = this.getOFD();//draw the angles each delta-angles to obtain the OFD
			double[] statistics = this.getStatistics(resultsOFD);//Analyse the OFD

			String titleText = "Angle: " + String.format( "%.2f", statistics[2]);
			titleText        = titleText+" Peak Density: "+String.format( "%.2f", statistics[1]);
			this.caller.plot.addText( titleText,0,0);
			this.caller.plot.draw();
        	this.caller.plot.show();
		}
		//create the OFD from the FFT
		double[] getOFD(){
			ImagePlus our_img_plus           = this.caller.ourImage;
			ImageProcessor our_img_processor = our_img_plus.getProcessor();	
			int size =  (int)Math.round(180/this.degree_separation)+1;//number of OFD elements
			double[] sumatoriesRatios = new double[size] ;
			double[] degrees   		  = new double[size] ;
			int count = 0;
			boolean drawLineAtAngle = this.caller.readCheckBox(DRAW_ANGLE);
			//for each delta-theta calculate the sumation of the fouriers transformation
			for(int degree = 0; degree<= 180; degree += this.degree_separation){
				double sumatoryRatio   = this.sumateAndDrawLineAtDegree(our_img_plus,our_img_processor,degree,drawLineAtAngle);
				sumatoriesRatios[count] = sumatoryRatio;
				degrees[count]	        = (double)degree;
				count++;
			}
			our_img_plus.show();
        	our_img_plus.updateAndDraw();
        	sumatoriesRatios = this.smooth(sumatoriesRatios,5);
        	//Plot the OFD
        	this.caller.plot = new Plot("OFD", "x", "y", degrees ,sumatoriesRatios ); 
        	this.caller.plot.draw();
        	this.caller.plot.show();
        	return sumatoriesRatios;
		}

		double sumateAndDrawLineAtDegree(ImagePlus our_img_plus,ImageProcessor our_img_processor,int angle,boolean drawLineAtAngle){
			//For each angle line the pixels are added
			double tangent = Math.tan(Math.toRadians(angle));
			int width      = our_img_processor.getWidth();
			int height     = our_img_processor.getHeight();

			int middle_width   = Math.round(width/2); 
			int middle_height  = Math.round(height/2);
			double sumOfPixels = 0.0;
			int numOfElements  = 0;
			//when the tangent is infinite is because the line is at 90 degrees
			if(tangent > 1E15){ 
				//Sum and draw vertical line at 90 degrees
				for(int i = 0;i<= middle_height;i++){
					double sumToTheVal  = (double)our_img_processor.get(middle_width,i);
					sumOfPixels += sumToTheVal/255;
					if(drawLineAtAngle)
						our_img_processor.set(middle_width,i,128);
					numOfElements ++;
				}
			}else{
				int start,end; 

				if( tangent >= 0 ){
					start = middle_width;
					end   = width;
				}else{
					start = 0;
					end   = middle_width;
				}
				//The line is define by y =mx+b, m = tan
				for(double x = start; x < end; x += 0.1  ){
					int y = Math.round((float)middle_height - (float)( (x-middle_width) * tangent) );
					if (y > 0 && y<height){
						double sumToTheVal  = (double)our_img_processor.get((int)Math.round(x),y);
						sumOfPixels += sumToTheVal/255;
						if(drawLineAtAngle)
							our_img_processor.set((int)Math.round(x),y,128);
						numOfElements ++;
					}
				}
			}
			if(drawLineAtAngle){
				our_img_plus.updateAndDraw();
				our_img_plus.draw();
				our_img_plus.show();
			}
			return (double)sumOfPixels/numOfElements;
		}
		// Obtain the peak and peak density from the OFD
		double[] getStatistics(double[] anglesSumatories)
		{
			int length = anglesSumatories.length;
			double max = Double.NEGATIVE_INFINITY;
			double min = Double.POSITIVE_INFINITY;
			double middle;
			double stepSize = 180.0/length;
			int indexOfMax = -1;
			//Obtain minumum and maximum values
			for(int i = 0; i< length ; i++){
				double ourVal = anglesSumatories[i];
				if(ourVal < min)
					min = ourVal;
				if(ourVal > max){
					max        = ourVal;
					indexOfMax = i;	
				}
			}
			middle = ((max-min)/2) + min;//Parameter to obtain the middle half
			//calculate the areas below the curve an peak to peak-density analysis
			double totalAreaBelow             = this.trapz(anglesSumatories, stepSize);
			double[] areaBelowThePeakandWidth = this.getAreaBelowThePeakAndWidthRatio(anglesSumatories,indexOfMax,middle,0.05);
			double areaBelowThePeak      = areaBelowThePeakandWidth[0];
			double widthRatioInPixels    = areaBelowThePeakandWidth[1];
			double areasRatio  = areaBelowThePeak/totalAreaBelow;
			double peakDensity = areasRatio/widthRatioInPixels;
			double angleOri    = indexOfMax * stepSize;
			double[] statistics = {max, peakDensity, angleOri};

			IJ.write("Angle Orientation [degrees]: " +  String.format( "%.2f", angleOri )  );
			IJ.write("Peak Density is: " + String.format( "%.2f", peakDensity )  );
			
			return statistics;
		}
		//method to smooth the OFD.	
		double[] smooth(double[] arrayToSmooth,int sizeOfSmooth){
			int  len     = arrayToSmooth.length;
			int   offset = (int)Math.ceil(sizeOfSmooth/2);
			double[] res = new double[len];
			for(int i =offset; i<len-offset; i++){
				double value = 0.0;
				for(int j = 0;j<sizeOfSmooth ;j++){
					int index = i-offset+j;
					value += arrayToSmooth[index];
				}
				res[i] = value/sizeOfSmooth;
			}
			return res;
		}
		//Obtain the area below the peak using the trapz method
		double[] getAreaBelowThePeakAndWidthRatio(double[] values,int indexOfMax,double middle,double tolerance)
		{
			int rigthIndex,leftIndex;
			double[] rigthSideValues,leftSideValues;
			boolean haveRigthLimit,haveLeftLimit = true;

			leftSideValues  = Arrays.copyOfRange(values,0,indexOfMax);
			rigthSideValues = Arrays.copyOfRange(values,indexOfMax,values.length);
			rigthIndex = this.getArrayClosestToBegingIndex(rigthSideValues,middle, tolerance);
			leftIndex  = this.getArrayClosestToEndIndex(leftSideValues,middle, tolerance);

			haveRigthLimit = rigthIndex != -1;
			haveLeftLimit  = leftIndex  != -1;

			if(!haveRigthLimit)
				rigthIndex =0;
			
			if(!haveLeftLimit)
				leftIndex =0;
			
			double area      = 0.0;
			int widthInteger = rigthIndex + (indexOfMax -  leftIndex ); 

			area = trapz(Arrays.copyOfRange(values,leftIndex, indexOfMax+rigthIndex), 180.0/values.length );

			if(!haveRigthLimit){//find the rigth limit in the other side of the OFD
				rigthIndex = this.getArrayClosestToBegingIndex(leftSideValues,middle, tolerance);
				area +=  trapz(Arrays.copyOfRange(leftSideValues,0, rigthIndex), 180/rigthIndex ); 
				widthInteger += rigthIndex;
			}
			if(!haveLeftLimit){//find the left limit in the other side of the OFD
				leftIndex = this.getArrayClosestToEndIndex(rigthSideValues,middle, tolerance);
				double localWidth = (rigthSideValues.length-leftIndex);
				area +=  trapz(Arrays.copyOfRange(rigthSideValues,leftIndex, rigthSideValues.length), 180/localWidth) ;
				widthInteger += localWidth;
			}
			double widthRatio   =  widthInteger/((double)values.length);
			double[] toReturn   =  {area, widthRatio };
			return toReturn;
		}
		//Helper to find the half width at half maximum
		public int getArrayClosestToBegingIndex(double[] arr,double value,double tolerance) {
		        int k = -1;
		        for(int i=0;i<arr.length;i++){
		            if(Math.abs(arr[i]-value) < tolerance){
		                k=i;
		                break;
		            }
		        }
		    return k;
		}
		//Helper to find the  width at half maximum
		public int getArrayClosestToEndIndex(double[] arr,double value,double tolerance) {
		        int k = -1;
		        for(int i=arr.length-1; i>=0 ;i--){
		            if(Math.abs(arr[i]-value) < tolerance ){
		                k=i;
		                break;
		            }
		        }
		    return k;
		}
		//Calculate the area below a discrete curve
		double trapz(double[] values,double h ){
			double integral = 0;
			int length = values.length;
			for(int i = 0; i<length-1; i++ )
				integral += h * 0.5*(values[i]+ values[i+1]);
			return integral;
		}
		void storeOriginalFft(){//only to keep track of the original FFT to turn back if required	
			ImagePlus ourImagePlus = this.caller.ourImage;
			ImagePlus fftBackup       = NewImage.createByteImage("Original Fft", ourImagePlus.getWidth(), ourImagePlus.getHeight(),
                                                     1, NewImage.FILL_BLACK);
			ImageProcessor fftBackupProcessor = fftBackup.getProcessor();
			fftBackupProcessor.copyBits(ourImagePlus.getProcessor(),0,0,Blitter.COPY);
			this.caller.fftImage = fftBackup;
		}
		void restore_original_fft(){//Get back to original FFT
			ImagePlus ourImagePlus 			 = this.caller.ourImage;
			ImageProcessor ourImageProcessor = ourImagePlus.getProcessor();
			ourImageProcessor.copyBits(this.caller.fftImage.getProcessor(),0,0,Blitter.COPY);
			ourImagePlus.show();
        	ourImagePlus.updateAndDraw();
		}
		void erode(){//Erode process to filter imagess
			ImagePlus ourImagePlus           = this.caller.ourImage;
			ImageProcessor ourImageProcessor = ourImagePlus.getProcessor();
			ourImageProcessor.dilate();
			ourImagePlus.show();
        	ourImagePlus.updateAndDraw();
		}
		void dilate(){//Erode process to filter imagess
			ImagePlus ourImagePlus 			 = this.caller.ourImage;
			ImageProcessor ourImageProcessor = ourImagePlus.getProcessor();
			ourImageProcessor.erode();
			ourImagePlus.show();
        	ourImagePlus.updateAndDraw();
		}
		//FFT transform and rotation of the image
		ImagePlus getFftandRotate(ImagePlus imp, ImageProcessor ip){
			FFT fft_transform = new FFT();
			fft_transform.run("");

			int w = ip.getWidth();
            int h = ip.getHeight();

            ImagePlus         fft_image = WindowManager.getCurrentImage();
            ImageProcessor fft_image_ip = fft_image.getProcessor();

        	fft_image_ip.copyBits(fft_image_ip,0,0,Blitter.COPY);
        	this.convertImageToGray8(fft_image);

        	ImageProcessor rotatedFftIP = fft_image_ip.rotateRight();
        	ImagePlus fft_image_rotated = new ImagePlus("Rotated FFT",rotatedFftIP);

        	fft_image_rotated.show();
        	fft_image_rotated.updateAndDraw();
        	fft_image.close();
        	return fft_image_rotated;
		}
		ImagePlus threshold(ImagePlus image_plus){//BINARIZE THE IMAGE
			this.convertImageToGray8(image_plus);
			ImageProcessor binarized_ip = image_plus.getProcessor();
			int h = binarized_ip.getHeight();
			int w = binarized_ip.getWidth();
			boolean binaryThreshold = this.caller.readCheckBox(BINARY_THRESHOLD);
			 for(int i =0;i<h;i++){
			 	for(int j=0; j<w; j++){
			 		int px = ((binarized_ip.get(j,i) ));
			 		if(px >this.threshold_value){
			 			if(binaryThreshold)
			 				binarized_ip.set(j,i,255);
			 		}else{
			 			binarized_ip.set(j,i,0);
			 		}
			 	}
			 }
			image_plus.show();
        	image_plus.updateAndDraw();
        	return image_plus;
		}
		void convertImageToGray8(ImagePlus img){//Basic image conversion
			ImageConverter image_converter  = new ImageConverter(img);
			image_converter.convertToGray8();
		}
	}
}
