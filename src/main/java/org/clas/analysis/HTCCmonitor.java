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
        this.setAnalysisTabNames("Electrons", "Reconstructed Hits");
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
        H2F hi_nphe_phi = new H2F("hi_nphe_phi", "hi_nphe_phi", 100, 0.0, 60.0, 100, -180.0, 180.0);  
        hi_nphe_phi.setTitleX("N.PhE"); 
        hi_nphe_phi.setTitleY("#phi_HTCC (deg)");
        H2F hi_nphe_theta = new H2F("hi_nphe_theta", "hi_nphe_theta", 100, 0.0, 60.0, 100, 0.0, 45.0);  
        hi_nphe_theta.setTitleX("N.PhE"); 
        hi_nphe_theta.setTitleY("#theta_HTCC (deg)");
        H1F hi_nphe_ele = new H1F("hi_nphe_ele", "hi_nphe_ele", 60, 0.0, 60.0);   
        hi_nphe_ele.setTitleX("N.PhE"); 
        hi_nphe_ele.setTitleY("Counts"); 
        H1F hi_nphe_cut = new H1F("hi_nphe_cut", "hi_nphe_cut", 60, 0.0, 60.0);   
        hi_nphe_cut.setTitleX("N.PhE"); 
        hi_nphe_cut.setTitleY("Counts"); 
        hi_nphe_cut.setLineColor(2);
        H2F hi_nphe_EC  = new H2F("hi_nphe_EC", "hi_nphe_EC", 100, 0, 0.6, 60, 0, 60);  
        hi_nphe_EC.setTitleX("E/p"); 
        hi_nphe_EC.setTitleY("N. PhE");
        H2F hi_phi_HBT = new H2F("hi_phi_HBT", "hi_phi_HBT", 100, -180.0, 180.0, 100, -180.0, 180.0);  
        hi_phi_HBT.setTitleX("#phi_HBT (deg)"); 
        hi_phi_HBT.setTitleY("#phi_HTCC (deg)");
        H2F hi_theta_HBT = new H2F("hi_theta_HBT", "hi_theta_HBT", 100, 0.0, 45.0, 100, 0.0, 45.0);  
        hi_theta_HBT.setTitleX("#theta_HBT (deg)"); 
        hi_theta_HBT.setTitleY("#theta_HTCC (deg)");
           
        DataGroup dg_rec = new DataGroup(2,2);
        dg_rec.addDataSet(hi_nphe, 0);
        dg_rec.addDataSet(hi_nhits, 1);
        dg_rec.addDataSet(hi_nphe_phi, 2);
        dg_rec.addDataSet(hi_nphe_theta, 3);
        this.getDataGroup().add(dg_rec,1);           
        DataGroup dg_elec = new DataGroup(2,2);
        dg_elec.addDataSet(hi_nphe_ele, 0);
        dg_elec.addDataSet(hi_nphe_cut, 0);
        dg_elec.addDataSet(hi_nphe_EC, 1);
        dg_elec.addDataSet(hi_phi_HBT, 2);
        dg_elec.addDataSet(hi_theta_HBT, 3);
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
        this.getAnalysisCanvas().getCanvas("Electrons").draw(this.getDataGroup().getItem(2).getH1F("hi_nphe_ele"));
        this.getAnalysisCanvas().getCanvas("Electrons").draw(this.getDataGroup().getItem(2).getH1F("hi_nphe_cut"),"same");
        this.getAnalysisCanvas().getCanvas("Electrons").cd(1);
        this.getAnalysisCanvas().getCanvas("Electrons").draw(this.getDataGroup().getItem(2).getH2F("hi_nphe_EC"));
        this.getAnalysisCanvas().getCanvas("Electrons").cd(2);
        this.getAnalysisCanvas().getCanvas("Electrons").draw(this.getDataGroup().getItem(2).getH2F("hi_phi_HBT"));
        this.getAnalysisCanvas().getCanvas("Electrons").cd(3);
        this.getAnalysisCanvas().getCanvas("Electrons").draw(this.getDataGroup().getItem(2).getH2F("hi_theta_HBT"));
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
        DataBank recBankEB = event.getBank("RECHB::Particle");
        DataBank recDeteEB = event.getBank("RECHB::Detector");

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
                
                if(recBankEB!=null && recDeteEB!=null) {
                    int nrows = recBankEB.rows();
                    for(int part = 0; part < nrows; part++){
                        int pidCode = 0;
                        if(recBankEB.getByte("charge", loop)==-1 && recBankEB.getByte("status", loop)==1) { // electron candidate
                            pidCode = 11;
                            Particle recParticle = new Particle(
                                                    pidCode,
                                                    recBankEB.getFloat("px", loop),
                                                    recBankEB.getFloat("py", loop),
                                                    recBankEB.getFloat("pz", loop),
                                                    recBankEB.getFloat("vx", loop),
                                                    recBankEB.getFloat("vy", loop),
                                                    recBankEB.getFloat("vz", loop));
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
                            if(energy1>0 && energy4>0) {
                                double energy=(energy1+energy4+energy7)/0.24;
                                this.getDataGroup().getItem(2).getH1F("hi_nphe_ele").fill(nphe*1.0);
                                this.getDataGroup().getItem(2).getH2F("hi_phi_HBT").fill(Math.toDegrees(recParticle.phi()),Math.toDegrees(phi));
                                this.getDataGroup().getItem(2).getH2F("hi_theta_HBT").fill(Math.toDegrees(recParticle.theta()),Math.toDegrees(theta));
                                if(Math.abs(Math.toDegrees(recParticle.phi())-phi)<30.0) {
                                    this.getDataGroup().getItem(2).getH1F("hi_nphe_cut").fill(nphe*1.0);
                                    this.getDataGroup().getItem(2).getH2F("hi_nphe_EC").fill((energy1+energy4+energy7)/recParticle.p(),nphe*1.0);
                                }
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
