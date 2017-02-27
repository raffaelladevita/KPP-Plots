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
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.F1D;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;

/**
 *
 * @author devita
 */
public class ELmonitor extends AnalysisMonitor {
    

    public ELmonitor(String name) {
        super(name);
        this.setAnalysisTabNames("Monte Carlo", "Electron Tracks");
        this.init(false);
    }

    
    @Override
    public void createHistos() {
        // initialize canvas and create histograms
        this.setNumberOfEvents(0);
        H1F summary = new H1F("summary","summary",6,1,7);
        summary.setTitleX("sector");
        summary.setTitleY("DC hits");
        summary.setFillColor(33);
        DataGroup sum = new DataGroup(1,1);
        sum.addDataSet(summary, 0);
        this.setAnalysisSummary(sum);
        // negative tracks
        H1F hi_p_neg = new H1F("hi_p_neg", "hi_p_neg", 100, 0.0, 8.0);     
        hi_p_neg.setTitleX("p (GeV)");
        hi_p_neg.setTitleY("Counts");
        H1F hi_theta_neg = new H1F("hi_theta_neg", "hi_theta_neg", 100, 0.0, 40.0); 
        hi_theta_neg.setTitleX("#theta (deg)");
        hi_theta_neg.setTitleY("Counts");
        H1F hi_phi_neg = new H1F("hi_phi_neg", "hi_phi_neg", 100, -180.0, 180.0);   
        hi_phi_neg.setTitleX("#phi (deg)");
        H1F hi_vz_neg = new H1F("hi_vz_neg", "hi_vz_neg", 100, -15.0, 15.0);   
        hi_vz_neg.setTitleX("Vz (cm)");
        hi_vz_neg.setTitleY("Counts");
        H1F hi_vz_neg_cut = new H1F("hi_vz_neg_cut", "hi_vz_neg_cut", 100, -15.0, 15.0);   
        hi_vz_neg_cut.setTitleX("Vz (cm)");
        hi_vz_neg_cut.setTitleY("Counts");
        hi_vz_neg_cut.setLineColor(2);
        F1D f1_vz_neg = new F1D("f1_vz_neg","[amp]*gaus(x,[mean],[sigma])", -8.0, 8.0);
        f1_vz_neg.setParameter(0, 0);
        f1_vz_neg.setParameter(1, 0);
        f1_vz_neg.setParameter(2, 1.0);
        f1_vz_neg.setLineWidth(2);
        f1_vz_neg.setLineColor(1);
        f1_vz_neg.setOptStat("1111");
        H2F hi_theta_p_neg = new H2F("hi_theta_p_neg", "hi_theta_p_neg", 100, 0.0, 8.0, 100, 0.0, 40.0); 
        hi_theta_p_neg.setTitleX("p (GeV)");
        hi_theta_p_neg.setTitleY("#theta (deg)");        
        H2F hi_theta_phi_neg = new H2F("hi_theta_phi_neg", "hi_theta_phi_neg", 100, -180.0, 180.0, 100, 0.0, 40.0); 
        hi_theta_phi_neg.setTitleX("#phi (deg)");
        hi_theta_phi_neg.setTitleY("#theta (deg)");        
        DataGroup dg_neg = new DataGroup(3,2);
        dg_neg.addDataSet(hi_p_neg, 0);
        dg_neg.addDataSet(hi_theta_neg, 1);
        dg_neg.addDataSet(hi_phi_neg, 2);
        dg_neg.addDataSet(hi_vz_neg, 3);
        dg_neg.addDataSet(hi_vz_neg_cut, 3);
        dg_neg.addDataSet(f1_vz_neg, 3);
        dg_neg.addDataSet(hi_theta_p_neg, 4);
        dg_neg.addDataSet(hi_theta_phi_neg, 5);
        this.getDataGroup().add(dg_neg, 1);        
        // mc comparison
        H1F hi_dp_neg = new H1F("hi_dp_neg", "hi_dp_neg", 100, -0.8, 0.8); 
        hi_dp_neg.setTitleX("#Delta P (GeV)");
        hi_dp_neg.setTitleY("Counts");
        hi_dp_neg.setTitle("Negative Tracks");
        H1F hi_dtheta_neg = new H1F("hi_dtheta_neg","hi_dtheta_neg", 100, -4.0, 4.0); 
        hi_dtheta_neg.setTitleX("#Delta #theta (deg)");
        hi_dtheta_neg.setTitleY("Counts");
        hi_dtheta_neg.setTitle("Negative Tracks");
        H1F hi_dphi_neg = new H1F("hi_dphi_neg", "hi_dphi_neg", 100, -8.0, 8.0); 
        hi_dphi_neg.setTitleX("#Delta #phi (deg)");
        hi_dphi_neg.setTitleY("Counts");
        hi_dphi_neg.setTitle("Negative Tracks");
        H1F hi_dvz_neg = new H1F("hi_dvz_neg", "hi_dvz_neg", 100, -20.0, 20.0);   
        hi_dvz_neg.setTitleX("#Delta Vz (cm)");
        hi_dvz_neg.setTitleY("Counts");
        hi_dvz_neg.setTitle("Negative Tracks");
        DataGroup mc = new DataGroup(2,2);
        mc.addDataSet(hi_dp_neg, 0);
        mc.addDataSet(hi_dtheta_neg, 1);
        mc.addDataSet(hi_dphi_neg, 2);
        mc.addDataSet(hi_dvz_neg, 3);
        this.getDataGroup().add(mc, 3);
    }
        
    @Override
    public void plotHistos() {
        this.getAnalysisCanvas().getCanvas("Electron Tracks").divide(3,2);
        this.getAnalysisCanvas().getCanvas("Electron Tracks").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Electron Tracks").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").divide(2, 2);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").setGridY(false);        this.getAnalysisCanvas().getCanvas("Electron Tracks").cd(0);
        this.getAnalysisCanvas().getCanvas("Electron Tracks").draw(this.getDataGroup().getItem(1).getH1F("hi_p_neg"));
        this.getAnalysisCanvas().getCanvas("Electron Tracks").cd(1);
        this.getAnalysisCanvas().getCanvas("Electron Tracks").draw(this.getDataGroup().getItem(1).getH1F("hi_theta_neg"));
        this.getAnalysisCanvas().getCanvas("Electron Tracks").cd(2);
        this.getAnalysisCanvas().getCanvas("Electron Tracks").draw(this.getDataGroup().getItem(1).getH1F("hi_phi_neg"));
        this.getAnalysisCanvas().getCanvas("Electron Tracks").cd(3);
        this.getAnalysisCanvas().getCanvas("Electron Tracks").draw(this.getDataGroup().getItem(1).getH1F("hi_vz_neg"));
        this.getAnalysisCanvas().getCanvas("Electron Tracks").draw(this.getDataGroup().getItem(1).getH1F("hi_vz_neg_cut"),"same");
        this.getAnalysisCanvas().getCanvas("Electron Tracks").draw(this.getDataGroup().getItem(1).getF1D("f1_vz_neg"),"same");
        this.getAnalysisCanvas().getCanvas("Electron Tracks").cd(4);
        this.getAnalysisCanvas().getCanvas("Electron Tracks").draw(this.getDataGroup().getItem(1).getH2F("hi_theta_p_neg"));
        this.getAnalysisCanvas().getCanvas("Electron Tracks").cd(5);
        this.getAnalysisCanvas().getCanvas("Electron Tracks").draw(this.getDataGroup().getItem(1).getH2F("hi_theta_phi_neg"));
        this.getAnalysisCanvas().getCanvas("Electron Tracks").update();
        this.getAnalysisCanvas().getCanvas("Monte Carlo").cd(0);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(this.getDataGroup().getItem(3).getH1F("hi_dp_neg"));
        this.getAnalysisCanvas().getCanvas("Monte Carlo").cd(1);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(this.getDataGroup().getItem(3).getH1F("hi_dtheta_neg"));
        this.getAnalysisCanvas().getCanvas("Monte Carlo").cd(2);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(this.getDataGroup().getItem(3).getH1F("hi_dphi_neg"));
        this.getAnalysisCanvas().getCanvas("Monte Carlo").cd(3);
        this.getAnalysisCanvas().getCanvas("Monte Carlo").draw(this.getDataGroup().getItem(3).getH1F("hi_dvz_neg"));
        this.getAnalysisCanvas().getCanvas("Electron Tracks").update();
        this.getAnalysisCanvas().getCanvas("Monte Carlo").update();
    }
    
    @Override
    public void processEvent(DataEvent event) {
        // process event info and save into data group
        Particle partGenNeg = null;
        Particle partRecNeg = null;
        if(event.hasBank("REC::Particle")==true){
            DataBank  bank = event.getBank("REC::Particle");
            int rows = bank.rows();
            for(int loop = 0; loop < rows; loop++){
                if(bank.getInt("pid", loop)==11 && partRecNeg==null) {
                    partRecNeg = new Particle(
                                          bank.getInt("pid", loop),
                                          bank.getFloat("px", loop),
                                          bank.getFloat("py", loop),
                                          bank.getFloat("pz", loop),
                                          bank.getFloat("vx", loop),
                                          bank.getFloat("vy", loop),
                                          bank.getFloat("vz", loop));
                    this.getDataGroup().getItem(1).getH1F("hi_p_neg").fill(partRecNeg.p());
                    this.getDataGroup().getItem(1).getH1F("hi_theta_neg").fill(Math.toDegrees(partRecNeg.theta()));
                    this.getDataGroup().getItem(1).getH1F("hi_phi_neg").fill(Math.toDegrees(partRecNeg.phi()));
                    this.getDataGroup().getItem(1).getH1F("hi_vz_neg").fill(partRecNeg.vz());
                    if(partRecNeg.p()>2.&& Math.toDegrees(partRecNeg.theta())>10.)this.getDataGroup().getItem(1).getH1F("hi_vz_neg_cut").fill(partRecNeg.vz());
                    this.getDataGroup().getItem(1).getH2F("hi_theta_p_neg").fill(partRecNeg.p(),Math.toDegrees(partRecNeg.theta()));
                    this.getDataGroup().getItem(1).getH2F("hi_theta_phi_neg").fill(Math.toDegrees(partRecNeg.phi()),Math.toDegrees(partRecNeg.theta()));                    
                }
            }
        }
        if(event.hasBank("MC::Particle")==true){
            DataBank genBank = event.getBank("MC::Particle");
            int nrows = genBank.rows();
            for(int loop = 0; loop < nrows; loop++) { 
                if(genBank.getInt("pid", loop)==11) {
                    partGenNeg = new Particle(
                                              genBank.getInt("pid", loop),
                                              genBank.getFloat("px", loop),
                                              genBank.getFloat("py", loop),
                                              genBank.getFloat("pz", loop),
                                              genBank.getFloat("vx", loop),
                                              genBank.getFloat("vy", loop),
                                              genBank.getFloat("vz", loop));
                    this.getDataGroup().getItem(3).getH1F("hi_dp_neg").fill(partRecNeg.p()-partGenNeg.p());
                    this.getDataGroup().getItem(3).getH1F("hi_dtheta_neg").fill(Math.toDegrees(partRecNeg.theta()-partGenNeg.theta()));
                    this.getDataGroup().getItem(3).getH1F("hi_dphi_neg").fill(Math.toDegrees(partRecNeg.phi()-partGenNeg.phi()));
                    this.getDataGroup().getItem(3).getH1F("hi_dvz_neg").fill(partRecNeg.vz()-partGenNeg.vz());
                }
            }
        }
    }

    @Override
    public void timerUpdate() {
//        System.out.println("Updating HBT");
        H1F hi_vz   = null;
        double mean  = 0;
        double sigma = 1.;
        double amp   = 100;
        //fitting negative tracks vertex
        hi_vz = this.getDataGroup().getItem(1).getH1F("hi_vz_neg_cut");
        mean  = hi_vz.getDataX(hi_vz.getMaximumBin());
        amp   = hi_vz.getBinContent(hi_vz.getMaximumBin());
        if(hi_vz.getEntries()<500) { // first fits 
            sigma = 1;
        }
        else {
            sigma = Math.abs(this.getDataGroup().getItem(1).getF1D("f1_vz_neg").getParameter(2));       
        }
        this.getDataGroup().getItem(1).getF1D("f1_vz_neg").setParameter(0, amp);
        this.getDataGroup().getItem(1).getF1D("f1_vz_neg").setParameter(1, mean);
        this.getDataGroup().getItem(1).getF1D("f1_vz_neg").setParameter(2, sigma);
        this.getDataGroup().getItem(1).getF1D("f1_vz_neg").setRange(mean-2.*sigma,mean+2.*sigma);
//        System.out.println(this.getDataGroup().getItem(1).getF1D("f1_vz_neg").getParameter(0) + " " + 
//                           this.getDataGroup().getItem(1).getF1D("f1_vz_neg").getParameter(1) + " " + 
//                           this.getDataGroup().getItem(1).getF1D("f1_vz_neg").getParameter(2));
        DataFitter.fit(this.getDataGroup().getItem(1).getF1D("f1_vz_neg"), hi_vz, "Q"); //No options uses error for sigma
    }
 

}
