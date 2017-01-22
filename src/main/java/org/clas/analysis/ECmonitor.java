/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.analysis;

import java.awt.BorderLayout;
import javax.swing.JSplitPane;
import org.clas.viewer.AnalysisMonitor;
import org.jlab.clas.physics.Particle;
import org.jlab.detector.base.DetectorType;
import org.jlab.detector.view.DetectorShape2D;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.graphics.EmbeddedCanvasTabbed;
import org.jlab.groot.group.DataGroup;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.rec.dc.Constants;
import org.jlab.rec.dc.GeometryLoader;

/**
 *
 * @author devita
 */
public class ECmonitor extends AnalysisMonitor {
    

    public ECmonitor(String name) {
        super(name);
        EmbeddedCanvasTabbed canvas = new EmbeddedCanvasTabbed("Electrons");
        canvas.addCanvas("Pizeros");
        this.setAnalysisCanvas(canvas);
        this.init();
    }

    
    @Override
    public void createHistos() {
        // initialize canvas and create histograms
        this.setNumberOfEvents(0);
        this.getAnalysisCanvas().getCanvas("Electrons").divide(2, 2);
        this.getAnalysisCanvas().getCanvas("Electrons").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Electrons").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Pizeros").divide(2, 2);
        this.getAnalysisCanvas().getCanvas("Pizeros").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Pizeros").setGridY(false);
        H1F summary = new H1F("summary","summary",6,1,7);
        summary.setTitleX("sector");
        summary.setTitleY("DC hits");
        summary.setFillColor(33);
        DataGroup sum = new DataGroup(1,1);
        sum.addDataSet(summary, 0);
        this.setAnalysisSummary(sum);
        // electrons
        H1F hi_dE_EC = new H1F("hi_dE_EC", "hi_dE_EC", 100, -0.8, 0.8);   
        hi_dE_EC.setTitleX("E-p (GeV)");
        hi_dE_EC.setTitleY("Counts");
        H1F hi_sf_EC = new H1F("hi_sf_EC", "hi_sf_EC", 100,  0.0, 0.6);   
        hi_sf_EC.setTitleX("E/p");
        hi_sf_EC.setTitleY("Counts");
        H2F hi_Evsp_EC = new H2F("hi_Evsp_EC", "hi_Evsp_EC", 100, 0, 6, 100, 0, 6);  
        hi_Evsp_EC.setTitleX("p (GeV)"); 
        hi_Evsp_EC.setTitleY("E (GeV)");
        H2F hi_sfvsp_EC = new H2F("hi_sfvsp_EC", "hi_sfvsp_EC", 100, 0, 6, 100, 0, 0.6);  
        hi_sfvsp_EC.setTitleX("p (GeV)"); 
        hi_sfvsp_EC.setTitleY("E/p");       
        DataGroup dg_electron = new DataGroup(2,2);
        dg_electron.addDataSet(hi_dE_EC, 0);
        dg_electron.addDataSet(hi_Evsp_EC, 1);
        dg_electron.addDataSet(hi_sf_EC, 2);
        dg_electron.addDataSet(hi_sfvsp_EC, 3);
        this.getDataGroup().add(dg_electron, 1);
        // pizero
        H1F hi_pi0_mass = new H1F("hi_pi0_mass","hi_pi0_mass",100,10.0,250.0);         
        hi_pi0_mass.setTitleX("Pizero Invariant Mass (MeV)");
        hi_pi0_mass.setTitleY("Counts");
        H1F hi_en_asy = new H1F("hi_en_asy","hi_en_asy",50,-1.0,1.0);      
        hi_en_asy.setTitleX("X:(E1-E2)/(E1+E2)");
        hi_en_asy.setTitleY("Counts");
        H1F hi_angle = new H1F("hi_angle","hi_angle",50, 0., 20.); 
        hi_angle.setTitleX("Two Photon Opening Angle (deg)");
        hi_angle.setTitleY("Counts");
        H2F hi_angle_en = new H2F("hi_angle_en","hi_angle_en",50, 0., 20., 50, 0, 1.5); 
        hi_angle_en.setTitleX("Two Photon Opening Angle (deg)");
        hi_angle_en.setTitleY("E1*E1 (GeV2)");
        DataGroup dg_pion = new DataGroup(2,2);
        dg_pion.addDataSet(hi_pi0_mass, 0);
        dg_pion.addDataSet(hi_en_asy, 1);
        dg_pion.addDataSet(hi_angle, 2);
        dg_pion.addDataSet(hi_angle_en, 3);
        this.getDataGroup().add(dg_pion, 2);
        
        this.getAnalysisCanvas().getCanvas("Electrons").cd(0);
        this.getAnalysisCanvas().getCanvas("Electrons").draw(hi_dE_EC);
        this.getAnalysisCanvas().getCanvas("Electrons").cd(1);
        this.getAnalysisCanvas().getCanvas("Electrons").draw(hi_sf_EC);
        this.getAnalysisCanvas().getCanvas("Electrons").cd(2);
        this.getAnalysisCanvas().getCanvas("Electrons").draw(hi_Evsp_EC);
        this.getAnalysisCanvas().getCanvas("Electrons").cd(3);
        this.getAnalysisCanvas().getCanvas("Electrons").draw(hi_sfvsp_EC);
        
        this.getAnalysisCanvas().getCanvas("Pizeros").cd(0);
        this.getAnalysisCanvas().getCanvas("Pizeros").draw(hi_pi0_mass);
        this.getAnalysisCanvas().getCanvas("Pizeros").cd(1);
        this.getAnalysisCanvas().getCanvas("Pizeros").draw(hi_en_asy);
        this.getAnalysisCanvas().getCanvas("Pizeros").cd(2);
        this.getAnalysisCanvas().getCanvas("Pizeros").draw(hi_angle);
        this.getAnalysisCanvas().getCanvas("Pizeros").cd(3);
        this.getAnalysisCanvas().getCanvas("Pizeros").draw(hi_angle_en);
        
        this.getAnalysisCanvas().getCanvas("Electrons").update();
        this.getAnalysisCanvas().getCanvas("Pizeros").update();

    }

    public void drawDetector() {
        // Load the Constants
        Constants.Load(true, true, 0);
        GeometryLoader.Load(10, "default");

    }

    @Override
    public void init() {
        this.getAnalysisPanel().setLayout(new BorderLayout());
        this.drawDetector();
        JSplitPane   splitPane = new JSplitPane();
        splitPane.setLeftComponent(this.getAnalysisView());
        splitPane.setRightComponent(this.getAnalysisCanvas());
        this.getAnalysisPanel().add(this.getAnalysisCanvas(),BorderLayout.CENTER);  
        this.createHistos();
    }
        
    @Override
    public void processEvent(DataEvent event) {
        // process event info and save into data group        
        Particle partRecEB = null;
        Particle partGamma1 = null;
        Particle partGamma2 = null;
        Particle partPi0    = null;
        // event builder
        DataBank recBankEB = event.getBank("REC::Particle");
        DataBank recDeteEB = event.getBank("REC::Detector");
        DataBank recBankTB = event.getBank("TimeBasedTrkg::TBTracks");
        if(recBankEB!=null && recDeteEB!=null) {
            int nrows = recBankEB.rows();
            for(int loop = 0; loop < nrows; loop++){
                    int pidCode = 0;
                    if(recBankEB.getByte("charge", loop)==-1) pidCode = 11;
                    else if(recBankEB.getByte("charge", loop)==1) pidCode = 211;
                    else pidCode = 22;
                    Particle recParticle = new Particle(
                                                pidCode,
                                                recBankEB.getFloat("px", loop),
                                                recBankEB.getFloat("py", loop),
                                                recBankEB.getFloat("pz", loop),
                                                recBankEB.getFloat("vx", loop),
                                                recBankEB.getFloat("vy", loop),
                                                recBankEB.getFloat("vz", loop));
                    if(recBankEB.getInt("pid", loop)==22) {
                        if(partGamma1==null)     partGamma1=recParticle;
                        else if(partGamma2==null) partGamma2=recParticle;
                    }
                    if(partRecEB==null && recBankEB.getByte("charge", loop)==-1 && recBankEB.getByte("status", loop)==1) {
                        recParticle.setProperty("sector",recBankTB.getByte("sector", loop)*1.0);
                        double energy1=0;
                        double energy4=0;
                        double energy7=0;
                        for(int j=0; j<recDeteEB.rows(); j++) {
                            if(recDeteEB.getShort("pindex",j)==loop && recDeteEB.getShort("detector",j)==16) {
                                if(energy1 == 0 && recDeteEB.getShort("layer",j) == 1) energy1 = recDeteEB.getFloat("energy",j);
                                if(energy4 == 0 && recDeteEB.getShort("layer",j) == 4) energy4 = recDeteEB.getFloat("energy",j);
                                if(energy7 == 0 && recDeteEB.getShort("layer",j) == 7) energy7 = recDeteEB.getFloat("energy",j);
                            }
                        }
                        recParticle.setProperty("energy1",energy1);
                        recParticle.setProperty("energy4",energy4);
                        recParticle.setProperty("energy7",energy7);
                        partRecEB=recParticle;
                        if(energy1>0 && energy4>0) {
                            double energy=(energy1+energy4+energy7)/0.24;
                            this.getDataGroup().getItem(1).getH1F("hi_dE_EC").fill(energy-partRecEB.p());
                            this.getDataGroup().getItem(1).getH1F("hi_sf_EC").fill((energy1+energy4+energy7)/partRecEB.p());
                            this.getDataGroup().getItem(1).getH2F("hi_Evsp_EC").fill(partRecEB.p(),energy);
                            this.getDataGroup().getItem(1).getH2F("hi_sfvsp_EC").fill(partRecEB.p(),(energy1+energy4+energy7)/partRecEB.p());
                        }
                    }
                }
                if(partGamma1!=null && partGamma2!=null) {
                    partPi0 = partGamma1;
                    partPi0.combine(partGamma2, +1);
                    double   invmass = 1e3*Math.sqrt(partPi0.vector().mass2());
                    double   x       = (partGamma1.p()-partGamma2.p())/(partGamma1.p()+partGamma2.p());
                    double   angle   = Math.toDegrees(Math.acos(partGamma1.cosTheta(partGamma2)));
                    this.getDataGroup().getItem(2).getH1F("hi_pi0_mass").fill(invmass);  //Two-photon invariant mass
                    this.getDataGroup().getItem(2).getH1F("hi_en_asy").fill(x);
                    this.getDataGroup().getItem(2).getH1F("hi_angle").fill(angle);
                    this.getDataGroup().getItem(2).getH2F("hi_angle_en").fill(angle,partGamma1.p()*partGamma2.p());
                }
        }
   }
   
   
    @Override
    public void resetEventListener() {
        System.out.println("Resetting EC histogram");
        this.createHistos();
    }

    @Override
    public void timerUpdate() {
//        System.out.println("Updating DC");
   }

    @Override
    public void setCanvasUpdate(int time) {
        this.getAnalysisCanvas().getCanvas("Electrons").initTimer(time);
        this.getAnalysisCanvas().getCanvas("Electrons").initTimer(time);
    }
 

}
