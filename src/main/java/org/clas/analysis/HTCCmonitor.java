/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.analysis;

import org.clas.viewer.AnalysisMonitor;
import org.jlab.clas.physics.Particle;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.group.DataGroup;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;

/**
 *
 * @author devita
 */
public class HTCCmonitor  extends AnalysisMonitor {
        
    
    public HTCCmonitor(String name) {
        super(name);
        this.setAnalysisTabNames("Reconstructed Hits", "Electrons");
        this.init(false);
    }

    @Override
    public void createHistos() {
        // initialize canvas and create histograms
        this.setNumberOfEvents(0);
        H1F summary = new H1F("summary","summary",6,1,7);
        summary.setTitleX("sector");
        summary.setTitleY("HTCC hits");
        summary.setFillColor(36);
        DataGroup sum = new DataGroup(1,1);
        sum.addDataSet(summary, 0);
        this.setAnalysisSummary(sum);
        H1F hi_nphe     = new H1F("hi_nphe", "hi_nphe", 60, 0.0, 60.0);   
        hi_nphe.setTitleX("N.PhE");
        hi_nphe.setTitleY("Counts");
        H1F hi_nhits     = new H1F("hi_nhits", "hi_nhits", 10, 0.0, 10.0);   
        hi_nhits.setTitleX("N.Hits in Rec. Clusters");
        hi_nhits.setTitleY("Counts");
        H2F hi_nphe_phi = new H2F("hi_nphe_phi", "hi_nphe_phi", 100, 0.0, 60.0, 24, -180.0, 180.0);  
        hi_nphe_phi.setTitleX("N.PhE"); 
        hi_nphe_phi.setTitleY("#phi_HTCC (deg)");
        H2F hi_nphe_theta = new H2F("hi_nphe_theta", "hi_nphe_theta", 100, 0.0, 60.0, 100, 0.0, 45.0);  
        hi_nphe_theta.setTitleX("N.PhE"); 
        hi_nphe_theta.setTitleY("#theta_HTCC (deg)");
        H1F hi_nphe_all = new H1F("hi_nphe_all", "hi_nphe_all", 60, 0.0, 60.0);   
        hi_nphe_all.setTitleX("N.PhE"); 
        hi_nphe_all.setTitleY("Counts"); 
        H1F hi_nphe_ele = new H1F("hi_nphe_ele", "hi_nphe_ele", 60, 0.0, 60.0);   
        hi_nphe_ele.setTitleX("N.PhE"); 
        hi_nphe_ele.setTitleY("Counts"); 
        hi_nphe_ele.setLineColor(2);
        H1F hi_time  = new H1F("hi_time", "hi_time", 100, -50., 50.);  
        hi_time.setTitleX("dT (ns)"); 
        hi_time.setTitleY("Counts");
        H2F hi_phi_TBT = new H2F("hi_phi_TBT", "hi_phi_TBT", 100, -180.0, 180.0, 100, -180.0, 180.0);  
        hi_phi_TBT.setTitleX("#phi_TBT (deg)"); 
        hi_phi_TBT.setTitleY("#phi_HTCC (deg)");
        H2F hi_theta_TBT = new H2F("hi_theta_TBT", "hi_theta_TBT", 100, 5.0, 45.0, 20, 5.0, 45.0);  
        hi_theta_TBT.setTitleX("#theta_TBT (deg)"); 
        hi_theta_TBT.setTitleY("#theta_HTCC (deg)");
           
        DataGroup dg_rec = new DataGroup(2,2);
        dg_rec.addDataSet(hi_nphe, 0);
        dg_rec.addDataSet(hi_nhits, 1);
        dg_rec.addDataSet(hi_nphe_phi, 2);
        dg_rec.addDataSet(hi_nphe_theta, 3);
        this.getDataGroup().add(dg_rec,1);           
        DataGroup dg_elec = new DataGroup(2,2);
        dg_elec.addDataSet(hi_nphe_all, 0);
        dg_elec.addDataSet(hi_nphe_ele, 0);
        dg_elec.addDataSet(hi_time, 1);
        dg_elec.addDataSet(hi_phi_TBT, 2);
        dg_elec.addDataSet(hi_theta_TBT, 3);
        this.getDataGroup().add(dg_elec,2);
    }
       
    @Override
    public void plotHistos() {
        this.getAnalysisCanvas().getCanvas("Electrons").divide(2, 2);
        this.getAnalysisCanvas().getCanvas("Electrons").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Electrons").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Reconstructed Hits").divide(2, 2);
        this.getAnalysisCanvas().getCanvas("Reconstructed Hits").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Reconstructed Hits").setGridY(false);
         
        // plotting histos
        this.getAnalysisCanvas().getCanvas("Electrons").cd(0);
        this.getAnalysisCanvas().getCanvas("Electrons").draw(this.getDataGroup().getItem(2).getH1F("hi_nphe_all"));
        this.getAnalysisCanvas().getCanvas("Electrons").draw(this.getDataGroup().getItem(2).getH1F("hi_nphe_ele"),"same");
        this.getAnalysisCanvas().getCanvas("Electrons").cd(1);
        this.getAnalysisCanvas().getCanvas("Electrons").draw(this.getDataGroup().getItem(2).getH1F("hi_time"));
        this.getAnalysisCanvas().getCanvas("Electrons").cd(2);
        this.getAnalysisCanvas().getCanvas("Electrons").draw(this.getDataGroup().getItem(2).getH2F("hi_phi_TBT"));
        this.getAnalysisCanvas().getCanvas("Electrons").cd(3);
        this.getAnalysisCanvas().getCanvas("Electrons").draw(this.getDataGroup().getItem(2).getH2F("hi_theta_TBT"));
        this.getAnalysisCanvas().getCanvas("Electrons").update();
        this.getAnalysisCanvas().getCanvas("Reconstructed Hits").cd(0);
        this.getAnalysisCanvas().getCanvas("Reconstructed Hits").draw(this.getDataGroup().getItem(1).getH1F("hi_nphe"));
        this.getAnalysisCanvas().getCanvas("Reconstructed Hits").cd(1);
        this.getAnalysisCanvas().getCanvas("Reconstructed Hits").draw(this.getDataGroup().getItem(1).getH1F("hi_nhits"));
        this.getAnalysisCanvas().getCanvas("Reconstructed Hits").cd(2);
        this.getAnalysisCanvas().getCanvas("Reconstructed Hits").draw(this.getDataGroup().getItem(1).getH2F("hi_nphe_phi"));
        this.getAnalysisCanvas().getCanvas("Reconstructed Hits").cd(3);
        this.getAnalysisCanvas().getCanvas("Reconstructed Hits").draw(this.getDataGroup().getItem(1).getH2F("hi_nphe_theta"));
        this.getAnalysisCanvas().getCanvas("Reconstructed Hits").update();
    }
    
    @Override
    public void processEvent(DataEvent event) {
        // process event info and save into data group
        DataBank recBankEB = null;
        DataBank recEvenEB = null;
        DataBank recDeteEB = null;
        if(event.hasBank("REC::Particle"))  recBankEB = event.getBank("REC::Particle");
        if(event.hasBank("REC::Event"))     recEvenEB = event.getBank("REC::Event");
        if(event.hasBank("REC::Cherenkov")) recDeteEB = event.getBank("REC::Cherenkov");
        if(event.hasBank("HTCC::rec")==true){
	    DataBank bank = event.getBank("HTCC::rec");
	    int rows = bank.rows();
	    for(int loop = 0; loop < rows; loop++){
                int  nphe   = bank.getInt("nphe", loop);
                int  nhits  = bank.getInt("nhits", loop);
                float phi   = bank.getFloat("phi", loop);
                float theta = bank.getFloat("theta", loop);
                this.getDataGroup().getItem(1).getH1F("hi_nphe").fill(nphe*1.0);
                this.getDataGroup().getItem(1).getH1F("hi_nhits").fill(nhits*1.0);
                this.getDataGroup().getItem(1).getH2F("hi_nphe_phi").fill(nphe*1.0,Math.toDegrees(phi));
                this.getDataGroup().getItem(1).getH2F("hi_nphe_theta").fill(nphe*1.0,Math.toDegrees(theta));
                
                if(recBankEB!=null && recDeteEB!=null && recEvenEB!=null) {
                    double startTime = recEvenEB.getFloat("STTime", 0);
                    int nrows = recBankEB.rows();
                    for(int part = 0; part < nrows; part++){
                        int pidCode = 0;
                        if(recBankEB.getInt("pid", loop)!=0) pidCode = recBankEB.getInt("pid", loop);
                        else if(recBankEB.getByte("charge", loop)==-1) pidCode = -211;
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
                        for(int j=0; j<recDeteEB.rows(); j++) {
                            if(recDeteEB.getShort("pindex",j)==loop && recDeteEB.getByte("detector",j)==15/*6*/) {
                                this.getDataGroup().getItem(2).getH1F("hi_nphe_all").fill(recDeteEB.getShort("nphe", j)*1.0);
                                if(recParticle.pid()==11) {
                                    this.getDataGroup().getItem(2).getH1F("hi_nphe_ele").fill(recDeteEB.getShort("nphe", j)*1.0);
                                    this.getDataGroup().getItem(2).getH1F("hi_time").fill(recDeteEB.getFloat("time", j)-recDeteEB.getFloat("path", j)/29.97-startTime);
                                }
                                this.getDataGroup().getItem(2).getH2F("hi_phi_TBT").fill(Math.toDegrees(recParticle.phi()),Math.toDegrees(recDeteEB.getFloat("phi", j)));
                                this.getDataGroup().getItem(2).getH2F("hi_theta_TBT").fill(Math.toDegrees(recParticle.theta()),Math.toDegrees(recDeteEB.getFloat("theta", j)));
                            }
                        }
                    }
                }
            }
    	}       
    }

    @Override
    public void timerUpdate() {

    }


}
