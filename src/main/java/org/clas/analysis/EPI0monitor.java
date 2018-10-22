/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.analysis;

import org.clas.viewer.AnalysisMonitor;
import org.jlab.clas.physics.LorentzVector;
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
public class EPI0monitor extends AnalysisMonitor {
    
    private double ebeam = 2.122193;
    
    public EPI0monitor(String name) {
        super(name);
        this.setAnalysisTabNames("Phi", "Electron", "MissingMass", "General");
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
        // General
        H2F hi_rec_q2w = new H2F("hi_rec_q2w","hi_rec_q2w",100, 0.6, ebeam*0.6+0.6, 100, 0., 1.); 
        hi_rec_q2w.setTitleX("W (GeV)");
        hi_rec_q2w.setTitleY("Q2 (GeV2)");
        H1F hi_rec_w = new H1F("hi_rec_w","hi_rec_w",250, 0.6, ebeam*0.6+0.6); 
        hi_rec_w.setTitleX("W (GeV)");
        hi_rec_w.setTitleY("Counts");
//        hi_rec_w.setOptStat("1110"); 
        F1D f1_w = new F1D("f1_w", "[amp]*gaus(x,[mean],[sigma])", 0.8, 1.2);
        f1_w.setParameter(0, 0);
        f1_w.setParameter(1, 1);
        f1_w.setParameter(2, 0.2);
        f1_w.setLineWidth(2);
        f1_w.setLineColor(2);
        f1_w.setOptStat("1111");
        H2F hi_rec_w_phi = new H2F("hi_rec_w_phi","hi_rec_w_phi",100, -180., 180., 250, 0.6, ebeam*0.6+0.6); 
        hi_rec_w_phi.setTitleX("#phi (deg)");
        hi_rec_w_phi.setTitleY("W (GeV)");
        H2F hi_rec_el = new H2F("hi_rec_el","hi_rec_el",100, 0.5, ebeam+0.5, 100, 0., 35.);
        hi_rec_el.setTitleX("p (GeV)");
        hi_rec_el.setTitleY("#theta (deg)");
        hi_rec_el.setTitle("Electron");
        H1F hi_rec_mm = new H1F("hi_rec_mm","hi_rec_mm",100, -0.5, 1);
        hi_rec_mm.setTitleX("Mx2 (GeV)");
        hi_rec_mm.setTitleY("Counts");
        hi_rec_mm.setTitle("MissingMass");       
        F1D f1_el = new F1D("f1_el", "2*(180/3.14)*atan(sqrt(0.93832*([e0]-x)/2/[e0]/x))", ebeam*0.75, ebeam*0.99);
        f1_el.setParameter(0, ebeam);
        H2F hi_w_mm = new H2F("hi_w_mm", "hi_w_mm", 100, 0.6, ebeam*0.6+0.6, 100,  -0.5, 1);
        hi_w_mm.setTitleX("W (GeV)");
        hi_w_mm.setTitleY("Mx2 (GeV)");
//        hi_dphi.setOptStat("1110");
        DataGroup dg_general = new DataGroup(2,3);
        dg_general.addDataSet(hi_rec_q2w,   0);
        dg_general.addDataSet(hi_rec_w,     1);
        dg_general.addDataSet(f1_w,         1);
        dg_general.addDataSet(hi_rec_w_phi, 2);
        dg_general.addDataSet(hi_rec_el,    3);
        dg_general.addDataSet(f1_el,        3);
        dg_general.addDataSet(hi_rec_mm,    4);  
        dg_general.addDataSet(hi_w_mm,      5);  
        this.getDataGroup().add(dg_general, 1);
        // Electron
        DataGroup dg_electron = new DataGroup(2,3);
        for(int sector=1; sector <= 6; sector++) {
            H1F hi_rec_w_sec = new H1F("hi_rec_w_" + sector, "hi_rec_w_" + sector, 250, 0.6, ebeam*0.6+0.6);  
            hi_rec_w_sec.setTitleX("W (GeV)");
            hi_rec_w_sec.setTitleY("Counts");
            hi_rec_w_sec.setTitle("Sector " + sector);
            F1D f1_w_sec = new F1D("f1_w_" + sector, "[amp]*gaus(x,[mean],[sigma])", 0.8, 1.2);
            f1_w_sec.setParameter(0, 0);
            f1_w_sec.setParameter(1, 1);
            f1_w_sec.setParameter(2, 0.2);
            f1_w_sec.setLineWidth(2);
            f1_w_sec.setLineColor(2);
            f1_w_sec.setOptStat("1111");
            dg_electron.addDataSet(hi_rec_w_sec, sector-1);
            dg_electron.addDataSet(f1_w_sec    , sector-1);
        }
        this.getDataGroup().add(dg_electron, 2);   
        // MissingMass
        DataGroup dg_proton = new DataGroup(2,3);
        for(int sector=1; sector <= 6; sector++) {
            H1F hi_rec_mm_sec = new H1F("hi_rec_mm_" + sector, "hi_rec_mm_" + sector, 100,  -0.5, 1.);
            hi_rec_mm_sec.setTitleX("Mx2 (GeV)");
            hi_rec_mm_sec.setTitleY("Counts");
            hi_rec_mm_sec.setTitle("Sector " + sector); 
            dg_proton.addDataSet(hi_rec_mm_sec, sector-1);  
        }
        this.getDataGroup().add(dg_proton, 3);   
        // Phi
        DataGroup dg_phi = new DataGroup(2,3);
        for(int sector=1; sector <= 6; sector++) {
            H1F hi_dphi_sec = new H1F("hi_dphi_" + sector, "hi_dphi_" + sector, 100, 140., 220.);  
            hi_dphi_sec.setTitleX("#Delta#phi (deg)");
            hi_dphi_sec.setTitleY("Counts");
            hi_dphi_sec.setTitle("Sector " + sector);
            F1D f1_dphi_sec = new F1D("f1_dphi_" + sector, "[amp]*gaus(x,[mean],[sigma])", 0.8, 1.2);
            f1_dphi_sec.setParameter(0, 0);
            f1_dphi_sec.setParameter(1, 1);
            f1_dphi_sec.setParameter(2, 0.2);
            f1_dphi_sec.setLineWidth(2);
            f1_dphi_sec.setLineColor(2);
            f1_dphi_sec.setOptStat("1111");
            dg_phi.addDataSet(hi_dphi_sec, sector-1);
            dg_phi.addDataSet(f1_dphi_sec, sector-1);
        }
        this.getDataGroup().add(dg_phi, 4);   
    }
        
    @Override
    public void plotHistos() {
        this.getAnalysisCanvas().getCanvas("General").divide(3,2);
        this.getAnalysisCanvas().getCanvas("General").setGridX(false);
        this.getAnalysisCanvas().getCanvas("General").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Electron").divide(3, 2);
        this.getAnalysisCanvas().getCanvas("Electron").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Electron").setGridY(false);
        this.getAnalysisCanvas().getCanvas("MissingMass").divide(3, 2);
        this.getAnalysisCanvas().getCanvas("MissingMass").setGridX(false);
        this.getAnalysisCanvas().getCanvas("MissingMass").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Phi").divide(3, 2);
        this.getAnalysisCanvas().getCanvas("Phi").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Phi").setGridY(false);
        this.getAnalysisCanvas().getCanvas("General").cd(0);
        this.getAnalysisCanvas().getCanvas("General").draw(this.getDataGroup().getItem(1).getH2F("hi_rec_q2w"));
        this.getAnalysisCanvas().getCanvas("General").getPad(0).getAxisZ().setLog(true);        
        this.getAnalysisCanvas().getCanvas("General").cd(1);
        this.getAnalysisCanvas().getCanvas("General").draw(this.getDataGroup().getItem(1).getH1F("hi_rec_w"));
        this.getAnalysisCanvas().getCanvas("General").draw(this.getDataGroup().getItem(1).getF1D("f1_w"),"same");
        this.getAnalysisCanvas().getCanvas("General").cd(2);
        this.getAnalysisCanvas().getCanvas("General").draw(this.getDataGroup().getItem(1).getH2F("hi_rec_w_phi"));
        this.getAnalysisCanvas().getCanvas("General").getPad(2).getAxisZ().setLog(true);        
        this.getAnalysisCanvas().getCanvas("General").cd(3);
        this.getAnalysisCanvas().getCanvas("General").draw(this.getDataGroup().getItem(1).getH2F("hi_rec_el"));
        this.getAnalysisCanvas().getCanvas("General").draw(this.getDataGroup().getItem(1).getF1D("f1_el"),"same");
        this.getAnalysisCanvas().getCanvas("General").getPad(3).getAxisZ().setLog(true);        
        this.getAnalysisCanvas().getCanvas("General").cd(4);
        this.getAnalysisCanvas().getCanvas("General").draw(this.getDataGroup().getItem(1).getH1F("hi_rec_mm"));
        this.getAnalysisCanvas().getCanvas("General").getPad(4).getAxisZ().setLog(true);        
        this.getAnalysisCanvas().getCanvas("General").cd(5);
        this.getAnalysisCanvas().getCanvas("General").draw(this.getDataGroup().getItem(1).getH2F("hi_w_mm"));
        this.getAnalysisCanvas().getCanvas("General").getPad(5).getAxisZ().setLog(true);        
        this.getAnalysisCanvas().getCanvas("General").update();
        for(int sector=1; sector <= 6; sector++) {
            this.getAnalysisCanvas().getCanvas("Electron").cd(sector-1);
            this.getAnalysisCanvas().getCanvas("Electron").draw(this.getDataGroup().getItem(2).getH1F("hi_rec_w_" + sector));
            this.getAnalysisCanvas().getCanvas("Electron").draw(this.getDataGroup().getItem(2).getF1D("f1_w_" + sector),"same");
        }
        this.getAnalysisCanvas().getCanvas("Electron").update();
         for(int sector=1; sector <= 6; sector++) {
            this.getAnalysisCanvas().getCanvas("MissingMass").cd(sector-1);
            this.getAnalysisCanvas().getCanvas("MissingMass").draw(this.getDataGroup().getItem(3).getH1F("hi_rec_mm_" + sector));
        }
        this.getAnalysisCanvas().getCanvas("MissingMass").update();        
         for(int sector=1; sector <= 6; sector++) {
            this.getAnalysisCanvas().getCanvas("Phi").cd(sector-1);
            this.getAnalysisCanvas().getCanvas("Phi").draw(this.getDataGroup().getItem(4).getH1F("hi_dphi_" + sector));
            this.getAnalysisCanvas().getCanvas("Phi").draw(this.getDataGroup().getItem(4).getF1D("f1_dphi_" + sector),"same");
        }
        this.getAnalysisCanvas().getCanvas("Phi").update();        
    }
    
    @Override
    public void processEvent(DataEvent event) {
        // process event info and save into data group
        Particle recEl = null;
        Particle recPr = null;
        LorentzVector virtualPhoton  = null;
        LorentzVector hadronSystem   = null;
        LorentzVector missingParticle  = null;
        if(event.hasBank("REC::Particle")==true && event.hasBank("REC::Track")){
            DataBank  bank  = event.getBank("REC::Particle");
            DataBank  track = event.getBank("REC::Track");
            int nCharged=0;
            int rows = bank.rows();
            for(int loop = 0; loop < rows; loop++){
                if(bank.getByte("charge", loop)!=0) nCharged++;
                if(bank.getInt("pid", loop)==11 && recEl==null && bank.getShort("status", loop)>=2000) {
                    recEl = new Particle(
                                          bank.getInt("pid", loop),
                                          bank.getFloat("px", loop),
                                          bank.getFloat("py", loop),
                                          bank.getFloat("pz", loop),
                                          bank.getFloat("vx", loop),
                                          bank.getFloat("vy", loop),
                                          bank.getFloat("vz", loop));
                    for(int j=0; j<track.rows(); j++) {
                        if(track.getShort("pindex", j)==loop) recEl.setProperty("sector", (double) track.getByte("sector", j));
                    }
                }
//                if(bank.getInt("charge", loop)==1 && recPr==null && bank.getShort("status", loop)<4000) {
                if(bank.getInt("pid", loop)==2212 && recPr==null && bank.getShort("status", loop)<4000) {
                    recPr = new Particle(
                                          2212,
                                          bank.getFloat("px", loop),
                                          bank.getFloat("py", loop),
                                          bank.getFloat("pz", loop),
                                          bank.getFloat("vx", loop),
                                          bank.getFloat("vy", loop),
                                          bank.getFloat("vz", loop));
                }
            }
            if(recEl != null && recPr != null && nCharged==2) {
//            System.out.println("Analyzed ");
                virtualPhoton = new LorentzVector(0., 0., ebeam, ebeam);
                virtualPhoton.sub(recEl.vector());
                hadronSystem = new LorentzVector(0., 0., ebeam, 0.9383+ebeam);
                hadronSystem.sub(recEl.vector());
                missingParticle = new LorentzVector(0., 0., ebeam, 0.9383+ebeam);
                missingParticle.sub(recEl.vector());
                missingParticle.sub(recPr.vector());
                int secEl = (int) recEl.getProperty("sector");
                if(Math.toDegrees(recEl.theta())>0){
                    this.getDataGroup().getItem(1).getH2F("hi_rec_q2w").fill(hadronSystem.mass(),-virtualPhoton.mass2());
                    this.getDataGroup().getItem(1).getH1F("hi_rec_w").fill(hadronSystem.mass());
                    this.getDataGroup().getItem(1).getH2F("hi_rec_w_phi").fill(Math.toDegrees(recEl.phi()), hadronSystem.mass());
                    this.getDataGroup().getItem(1).getH2F("hi_rec_el").fill(recEl.p(),Math.toDegrees(recEl.theta()));
                    this.getDataGroup().getItem(2).getH1F("hi_rec_w_" + secEl).fill(hadronSystem.mass());
                }
                this.getDataGroup().getItem(1).getH1F("hi_rec_mm").fill(missingParticle.mass2());
//                if(hadronSystem.mass()>1.1 && hadronSystem.mass()<1.4) {
                    this.getDataGroup().getItem(1).getH2F("hi_w_mm").fill(hadronSystem.mass(),missingParticle.mass2());
                    this.getDataGroup().getItem(3).getH1F("hi_rec_mm_" + secEl).fill(missingParticle.mass2());
                    this.getDataGroup().getItem(4).getH1F("hi_dphi_" + secEl).fill(Math.abs(Math.toDegrees(recPr.phi()-recEl.phi())));
//                }
            }

        }
    }

    @Override
    public void timerUpdate() {
        this.fitW(this.getDataGroup().getItem(1).getH1F("hi_rec_w"), this.getDataGroup().getItem(1).getF1D("f1_w"));
//        for(int sector=1; sector <= 6; sector++) {
//            this.fitW(this.getDataGroup().getItem(2).getH1F("hi_rec_w_" + sector), this.getDataGroup().getItem(2).getF1D("f1_w_" + sector));
//            this.fitPhi(this.getDataGroup().getItem(4).getH1F("hi_dphi_" + sector), this.getDataGroup().getItem(4).getF1D("f1_dphi_" + sector));
//        }
    }
    
    public void fitW(H1F hiw,F1D f1w) {

        double mean = hiw.getDataX(hiw.getMaximumBin());
        double amp = hiw.getBinContent(hiw.getMaximumBin());
        double sigma = 0.05;
        f1w.setParameter(0, amp);
        f1w.setParameter(1, mean);
        f1w.setParameter(2, sigma);
        double rmax = mean + 1. * Math.abs(sigma);
        double rmin = mean - 2. * Math.abs(sigma);
        f1w.setRange(rmin, rmax);
        DataFitter.fit(f1w, hiw, "Q"); //No options uses error for sigma 
        hiw.setFunction(null);
        mean = f1w.getParameter(1);
        sigma = f1w.getParameter(2);
        rmax = mean + 1. * Math.abs(sigma);
        rmin = mean - 2. * Math.abs(sigma);
        f1w.setRange(rmin, rmax);
//        System.out.println(mean + " " + sigma + " " + rmin + " " + rmax);
        DataFitter.fit(f1w, hiw, "Q"); //No options uses error for sigma 
        hiw.setFunction(null);
    }

    public void fitPhi(H1F hiw,F1D f1w) {

        double mean = hiw.getDataX(hiw.getMaximumBin());
        double amp = hiw.getBinContent(hiw.getMaximumBin());
        double sigma = 3;
        f1w.setParameter(0, amp);
        f1w.setParameter(1, mean);
        f1w.setParameter(2, sigma);
        double rmax = mean + 2. * Math.abs(sigma);
        double rmin = mean - 2. * Math.abs(sigma);
        f1w.setRange(rmin, rmax);
        DataFitter.fit(f1w, hiw, "Q"); //No options uses error for sigma 
        hiw.setFunction(null);
        mean = f1w.getParameter(1);
        sigma = f1w.getParameter(2);
        rmax = mean + 2. * Math.abs(sigma);
        rmin = mean - 2. * Math.abs(sigma);
        f1w.setRange(rmin, rmax);
//        System.out.println(mean + " " + sigma + " " + rmin + " " + rmax);
        DataFitter.fit(f1w, hiw, "Q"); //No options uses error for sigma 
        hiw.setFunction(null);
    }

}
