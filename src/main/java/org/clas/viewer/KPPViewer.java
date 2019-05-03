package org.clas.viewer;

import org.clas.analysis.HBTmonitor;
import java.awt.BorderLayout;
import java.awt.Font;
import java.awt.Image;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.TreeMap;
import javax.imageio.ImageIO;
import javax.swing.ImageIcon;
import javax.swing.JFileChooser;

import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextPane;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.text.SimpleAttributeSet;
import javax.swing.text.StyleConstants;
import javax.swing.text.StyledDocument;
import org.clas.analysis.CNDmonitor;
import org.clas.analysis.CTOFmonitor;
import org.clas.analysis.CVTmonitor;
import org.clas.analysis.EBHBmonitor;
import org.clas.analysis.EBmonitor;
import org.clas.analysis.ECmonitor;
import org.clas.analysis.ELASTICmonitor;
import org.clas.analysis.ELmonitor;
import org.clas.analysis.EPI0monitor;
import org.clas.analysis.EPIPLUSmonitor;
import org.clas.analysis.FMTmonitor;
import org.clas.analysis.FTOFmonitor;
import org.clas.analysis.FTmonitor;
import org.clas.analysis.HTCCmonitor;
import org.clas.analysis.KINEmonitor;
import org.clas.analysis.LTCCmonitor;
import org.clas.analysis.TBTmonitor;
import org.clas.analysis.TIMEmonitor;

import org.jlab.detector.decode.CodaEventDecoder;
import org.jlab.detector.decode.DetectorEventDecoder;
import org.jlab.detector.view.DetectorListener;
import org.jlab.detector.view.DetectorPane2D;
import org.jlab.detector.view.DetectorShape2D;
import org.jlab.groot.base.GStyle;
import org.jlab.groot.data.H2F;
import org.jlab.groot.data.TDirectory;
import org.jlab.groot.graphics.EmbeddedCanvasTabbed;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.base.DataEventType;
import org.jlab.io.task.DataSourceProcessorPane;
import org.jlab.io.task.IDataEventListener;

/**
 *
 * @author ziegler
 */
public class KPPViewer implements IDataEventListener, DetectorListener, ActionListener, ChangeListener {
    
    List<DetectorPane2D> AnalysisPanels 	= new ArrayList<DetectorPane2D>();
    JTabbedPane tabbedpane           		= null;
    JPanel mainPanel 				= null;
    JMenuBar menuBar                            = null;
    DataSourceProcessorPane processorPane 	= null;
    EmbeddedCanvasTabbed CLAS12Canvas           = null;

    
    CodaEventDecoder               decoder = new CodaEventDecoder();
    DetectorEventDecoder   detectorDecoder = new DetectorEventDecoder();
       
    
    TreeMap<String, List<H2F>>  histos = new TreeMap<String,List<H2F>>();
    
    private int canvasUpdateTime = 10000;
    private int analysisUpdateTime = 10000;
    private int runNumber  = 0;
    private String kppDir = "/Users/devita";
    
    // detector monitors
    AnalysisMonitor[] monitors = {
    		new HBTmonitor("HBT"),
    		new TBTmonitor("TBT"),
                new CVTmonitor("CVT"),
                new ELmonitor("ELECTRONS"),
                new ECmonitor("EC"),
                new FTOFmonitor("FTOF"),
        	new HTCCmonitor("HTCC"),
        	new LTCCmonitor("LTCC"),
                new CTOFmonitor("CTOF"),
                new CNDmonitor("CND"),
                new FTmonitor("FT"),
                new FMTmonitor("FMT"),
        	new EBHBmonitor("EBHB"),
        	new EBmonitor("EB"),
        	new KINEmonitor("KINEMATICS"),
           	new ELASTICmonitor("ELASTIC"),
        	new EPIPLUSmonitor("EPIPLUS"),
        	new EPI0monitor("EPI0"),
        	new TIMEmonitor("TIME")
    };
        
    public KPPViewer() {    	
        		
	// create menu bar
        menuBar = new JMenuBar();
        JMenuItem menuItem;
        JMenu file = new JMenu("File");
        file.getAccessibleContext().setAccessibleDescription("File options");
        menuItem = new JMenuItem("Open histograms file...");
        menuItem.getAccessibleContext().setAccessibleDescription("Open histograms file");
        menuItem.addActionListener(this);
        file.add(menuItem);
         menuItem = new JMenuItem("Print histograms to file...");
        menuItem.getAccessibleContext().setAccessibleDescription("Print histograms to file");
        menuItem.addActionListener(this);
        file.add(menuItem);
        menuItem = new JMenuItem("Save histograms to file...");
        menuItem.getAccessibleContext().setAccessibleDescription("Save histograms to file");
        menuItem.addActionListener(this);
        file.add(menuItem);
        menuBar.add(file);
        JMenu settings = new JMenu("Settings");
        settings.getAccessibleContext().setAccessibleDescription("Choose monitoring parameters");
        menuItem = new JMenuItem("Set GUI update interval...");
        menuItem.getAccessibleContext().setAccessibleDescription("Set GUI update interval");
        menuItem.addActionListener(this);
        settings.add(menuItem);
        menuItem = new JMenuItem("Set Style to default");
        menuItem.getAccessibleContext().setAccessibleDescription("Set GROOT style to default");
        menuItem.addActionListener(this);
        settings.add(menuItem);
        menuItem = new JMenuItem("Set Style for performance plots");
        menuItem.getAccessibleContext().setAccessibleDescription("Set GROOT style for performance plots");
        menuItem.addActionListener(this);
        settings.add(menuItem);
        menuBar.add(settings);
        
           
        // create main panel
        mainPanel = new JPanel();	
	mainPanel.setLayout(new BorderLayout());
        
      	tabbedpane 	= new JTabbedPane();

        processorPane = new DataSourceProcessorPane();
        processorPane.setUpdateRate(analysisUpdateTime);

        mainPanel.add(tabbedpane);
        mainPanel.add(processorPane,BorderLayout.PAGE_END);
        
    
        GStyle.getAxisAttributesX().setTitleFontSize(24);
        GStyle.getAxisAttributesX().setLabelFontSize(18);
        GStyle.getAxisAttributesY().setTitleFontSize(24);
        GStyle.getAxisAttributesY().setLabelFontSize(18);
        CLAS12Canvas    = new EmbeddedCanvasTabbed("Summaries");
        CLAS12Canvas.getCanvas("Summaries").divide(2,2);
        CLAS12Canvas.getCanvas("Summaries").setGridX(false);
        CLAS12Canvas.getCanvas("Summaries").setGridY(false);
        JPanel    CLAS12View = new JPanel(new BorderLayout());
        JSplitPane splitPanel = new JSplitPane();
        splitPanel.setLeftComponent(CLAS12View);
        splitPanel.setRightComponent(CLAS12Canvas);
        JTextPane clas12Text   = new JTextPane();
        clas12Text.setText("CLAS12\n KPP plots\n V3.0");
        clas12Text.setEditable(false);
        StyledDocument styledDoc = clas12Text.getStyledDocument();
        SimpleAttributeSet center = new SimpleAttributeSet();
        StyleConstants.setAlignment(center, StyleConstants.ALIGN_CENTER);
        styledDoc.setParagraphAttributes(0, styledDoc.getLength(), center, false);
        clas12Text.setBackground(CLAS12View.getBackground());
        clas12Text.setFont(new Font("Avenir",Font.PLAIN,20));
        JLabel clas12Design = this.getImage("https://www.jlab.org/Hall-B/clas12-web/sidebar/clas12-design.jpg",0.1);
        JLabel clas12Logo   = this.getImage("https://www.jlab.org/Hall-B/pubs-web/logo/CLAS-frame-low.jpg", 0.3);
//        CLAS12View.add(clas12Name,BorderLayout.PAGE_START);
        CLAS12View.add(clas12Design);
        CLAS12View.add(clas12Text,BorderLayout.PAGE_END);
 
        
        tabbedpane.add(splitPanel,"CLAS12");
        tabbedpane.addChangeListener(this);
       
        for(int k =0; k<this.monitors.length; k++) {
                this.tabbedpane.add(this.monitors[k].getAnalysisPanel(), this.monitors[k].getAnalysisName());
        	this.monitors[k].getAnalysisView().getView().addDetectorListener(this);
        }
        this.processorPane.addEventListener(this);
        
        this.setCanvasUpdate(canvasUpdateTime);
        this.plotSummaries();
    }
      
    public void actionPerformed(ActionEvent e) {
        System.out.println(e.getActionCommand());
        if(e.getActionCommand()=="Set GUI update interval...") {
            this.chooseUpdateInterval();
        }
        if(e.getActionCommand()=="Set Style to default") {
            this.setStyle(0);
        }
        if(e.getActionCommand()=="Set Style for performance plots") {
            this.setStyle(1);
        }
        if(e.getActionCommand()=="Open histograms file...") {
            String fileName = null;
            JFileChooser fc = new JFileChooser();
            fc.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
            File workingDirectory = new File(this.kppDir + "/kpp-histos");
            fc.setCurrentDirectory(workingDirectory);
            int option = fc.showOpenDialog(null);
            if (option == JFileChooser.APPROVE_OPTION) {
                fileName = fc.getSelectedFile().getAbsolutePath();            
            }
            if(fileName != null) this.loadHistosFromFile(fileName);
        }        
        if(e.getActionCommand()=="Print histograms to file...") {
            this.printHistosToFile();
        }
        if(e.getActionCommand()=="Save histograms to file...") {
            DateFormat df = new SimpleDateFormat("MM-dd-yyyy_hh.mm.ss_aa");
            String fileName = "clas12rec_run_" + this.runNumber + "_" + df.format(new Date()) + ".hipo";
            JFileChooser fc = new JFileChooser();
            File workingDirectory = new File(this.kppDir + "/kpp-histos");
            fc.setCurrentDirectory(workingDirectory);
            File file = new File(fileName);
            fc.setSelectedFile(file);
            int returnValue = fc.showSaveDialog(null);
            if (returnValue == JFileChooser.APPROVE_OPTION) {
               fileName = fc.getSelectedFile().getAbsolutePath();            
            }
            this.saveHistosToFile(fileName);
        }
    }

    public void chooseUpdateInterval() {
        String s = (String)JOptionPane.showInputDialog(
                    null,
                    "GUI update interval (ms)",
                    " ",
                    JOptionPane.PLAIN_MESSAGE,
                    null,
                    null,
                    "1000");
        if(s!=null){
            int time = 1000;
            try { 
                time= Integer.parseInt(s);
            } catch(NumberFormatException e) { 
                JOptionPane.showMessageDialog(null, "Value must be a positive integer!");
            }
            if(time>0) {
                this.setCanvasUpdate(time);
            }
            else {
                JOptionPane.showMessageDialog(null, "Value must be a positive integer!");
            }
        }
    }
        
    private JLabel getImage(String path,double scale) {
        JLabel label = null;
        Image image = null;
        try {
            URL url = new URL(path);
            image = ImageIO.read(url);
        } catch (IOException e) {
        	e.printStackTrace();
                System.out.println("Picture upload from " + path + " failed");
        }
        ImageIcon imageIcon = new ImageIcon(image);
        double width  = imageIcon.getIconWidth()*scale;
        double height = imageIcon.getIconHeight()*scale;
        imageIcon = new ImageIcon(image.getScaledInstance((int) width,(int) height, Image.SCALE_SMOOTH));
        label = new JLabel(imageIcon);
        return label;
    }
    
    public JPanel  getPanel(){
        return mainPanel;
    }

    private int getRunNumber(DataEvent event) {
        int rNum = this.runNumber;
        DataBank bank = event.getBank("RUN::config");
        if(bank!=null) {
            rNum = bank.getInt("run", 0);
        }
        return rNum;
    }
    
    @Override
    public void dataEventAction(DataEvent event) {
    	
       // EvioDataEvent decodedEvent = deco.DecodeEvent(event, decoder, table);
        //decodedEvent.show();
        		
	if(event!=null ){
//            event.show();
            if (event.getType() == DataEventType.EVENT_START) {
                this.runNumber = this.getRunNumber(event);
            }
            if(this.runNumber != this.getRunNumber(event)) {
//                this.saveToFile("mon12_histo_run_" + runNumber + ".hipo");
                this.runNumber = this.getRunNumber(event);
//                resetEventListener();
            }
            for(int k=0; k<this.monitors.length; k++) {
                this.monitors[k].dataEventAction(event);
            }      
            if (event.getType() == DataEventType.EVENT_STOP) {
                for(int k=0; k<this.monitors.length; k++) {
                    this.monitors[k].timerUpdate();
                }                
            }
	}
   }

    public void loadHistosFromFile(String fileName) {
        // TXT table summary FILE //
        System.out.println("Opening file: " + fileName);
        TDirectory dir = new TDirectory();
        dir.readFile(fileName);
        System.out.println(dir.getDirectoryList());
        dir.cd();
        dir.pwd();
        
        for(int k=0; k<this.monitors.length; k++) {
            this.monitors[k].readDataGroup(dir);
        }
        this.plotSummaries();
    }

    public void plotSummaries() {
        this.CLAS12Canvas.getCanvas("Summaries").cd(0);
        
        if(this.monitors[1].getDataGroup().getItem(1).getH1F("hi_vz_neg_cut")!=null) {
            this.CLAS12Canvas.getCanvas("Summaries").draw(this.monitors[1].getDataGroup().getItem(1).getH1F("hi_vz_neg"));
            this.CLAS12Canvas.getCanvas("Summaries").draw(this.monitors[1].getDataGroup().getItem(1).getF1D("f1_vz_neg"),"same");
        }
        this.CLAS12Canvas.getCanvas("Summaries").cd(1);
        if(this.monitors[4].getDataGroup().getItem(1).getH2F("hi_Evsp_EC")!=null) this.CLAS12Canvas.getCanvas("Summaries").draw(this.monitors[4].getDataGroup().getItem(1).getH2F("hi_Evsp_EC"));
        this.CLAS12Canvas.getCanvas("Summaries").cd(3);
        if(this.monitors[4].getDataGroup().getItem(1).getH2F("hi_sfvsp_EC")!=null) this.CLAS12Canvas.getCanvas("Summaries").draw(this.monitors[4].getDataGroup().getItem(1).getH2F("hi_sfvsp_EC"));
        this.CLAS12Canvas.getCanvas("Summaries").cd(2);
        if(this.monitors[4].getDataGroup().getItem(2).getH1F("hi_pi0_mass")!=null) {
                this.CLAS12Canvas.getCanvas("Summaries").draw(this.monitors[4].getDataGroup().getItem(2).getH1F("hi_pi0_mass"));
                this.CLAS12Canvas.getCanvas("Summaries").draw(this.monitors[4].getDataGroup().getItem(2).getF1D("fpi0"),"same");
        }    
    }
    
    public void printHistosToFile() {
        DateFormat df = new SimpleDateFormat("MM-dd-yyyy_hh.mm.ss_aa");
        String data = this.kppDir + "/clas12rec_run_" + this.runNumber + "_" + df.format(new Date());        
        File theDir = new File(data);
        // if the directory does not exist, create it
        if (!theDir.exists()) {
            boolean result = false;
            try{
                theDir.mkdir();
                result = true;
            } 
            catch(SecurityException se){
                //handle it
            }        
            if(result) {    
            System.out.println("Created directory: " + data);
            }
        }
        String fileName = data + "/clas12_canvas.png";
        System.out.println(fileName);
        this.CLAS12Canvas.getCanvas("Summaries").save(fileName);
        for(int k=0; k<this.monitors.length; k++) {
            this.monitors[k].printCanvas(data);
        }
    }

     @Override
    public void processShape(DetectorShape2D shape) {
        System.out.println("SHAPE SELECTED = " + shape.getDescriptor());
    }
    
    @Override
    public void resetEventListener() {
        for(int k=0; k<this.monitors.length; k++) {
            this.monitors[k].resetEventListener();
            this.monitors[k].timerUpdate();
        }      
        this.plotSummaries();
    }
    
   public void saveHistosToFile(String fileName) {
        // TXT table summary FILE //
        TDirectory dir = new TDirectory();
        for(int k=0; k<this.monitors.length; k++) {
            this.monitors[k].writeDataGroup(dir);
        }
        System.out.println("Saving histograms to file " + fileName);
        dir.writeFile(fileName);
    }
            
    public void setCanvasUpdate(int time) {
        System.out.println("Setting " + time + " ms update interval");
        this.canvasUpdateTime = time;
        this.CLAS12Canvas.getCanvas("Summaries").initTimer(time);
        this.CLAS12Canvas.getCanvas("Summaries").update();
        for(int k=0; k<this.monitors.length; k++) {
            this.monitors[k].setCanvasUpdate(time);
        }
    }

    public void setStyle(int mode) {
        for(int k=0; k<this.monitors.length; k++) {
            this.monitors[k].setStyle(mode);
            this.plotSummaries();
            this.monitors[k].plotHistos();
        }    
    }
    
    public void stateChanged(ChangeEvent e) {
        this.timerUpdate();
    }
    
    @Override
    public void timerUpdate() {
//        System.out.println("Time to update ...");
        for(int k=0; k<this.monitors.length; k++) {
            this.monitors[k].timerUpdate();
        }
    }

    public static void main(String[] args){
        JFrame frame = new JFrame("KPP-Plots");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        KPPViewer viewer = new KPPViewer();
        //frame.add(viewer.getPanel());
        frame.add(viewer.mainPanel);
        frame.setJMenuBar(viewer.menuBar);
        frame.setSize(1400, 800);
        frame.setVisible(true);
    }
   
}