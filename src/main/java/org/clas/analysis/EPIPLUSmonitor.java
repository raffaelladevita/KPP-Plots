/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.analysis;

import org.clas.viewer.AnalysisMonitor;
import org.jlab.clas.physics.LorentzVector;
import org.jlab.clas.physics.Particle;
import org.jlab.detector.calib.utils.ConstantsManager;
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
public class EPIPLUSmonitor extends AnalysisMonitor {
    
    private double ebeam = 10.6;
    
    public EPIPLUSmonitor(String name, ConstantsManager ccdb) {
        super(name,ccdb);
        this.setAnalysisTabNames("MissingMass", "General");
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
        H2F hi_rec_q2w = new H2F("hi_rec_q2w","hi_rec_q2w",100, 0.6, ebeam*0.4+0.6, 100, 0.5, 4.5); 
        hi_rec_q2w.setTitleX("W (GeV)");
        hi_rec_q2w.setTitleY("Q2 (GeV2)");
        H1F hi_rec_w = new H1F("hi_rec_w","hi_rec_w",250, 0.6, ebeam*0.4+0.6); 
        hi_rec_w.setTitleX("W (GeV)");
        hi_rec_w.setTitleY("Counts");
//        hi_rec_w.setOptStat("1110"); 
        H2F hi_rec_w_phi = new H2F("hi_rec_w_phi","hi_rec_w_phi",100, -180., 180., 250, 0.6, ebeam*0.4+0.6); 
        hi_rec_w_phi.setTitleX("#phi (deg)");
        hi_rec_w_phi.setTitleY("W (GeV)");
        H2F hi_rec_el = new H2F("hi_rec_el","hi_rec_el",100, 0.5, ebeam+0.5, 100, 0., 35.);
        hi_rec_el.setTitleX("p (GeV)");
        hi_rec_el.setTitleY("#theta (deg)");
        hi_rec_el.setTitle("Electron");
        H1F hi_rec_mm = new H1F("hi_rec_mm","hi_rec_mm",200, 0.6, ebeam*0.4+0.6);
        hi_rec_mm.setTitleX("Mx (GeV)");
        hi_rec_mm.setTitleY("Counts");
        hi_rec_mm.setTitle("MissingMass"); 
        F1D f1_mm = new F1D("f1_mm", "[amp]*gaus(x,[mean],[sigma])", 0.8, 1.2);
        f1_mm.setParameter(0, 0);
        f1_mm.setParameter(1, 1);
        f1_mm.setParameter(2, 0.2);
        f1_mm.setLineWidth(2);
        f1_mm.setLineColor(2);
        f1_mm.setOptStat("1111");        
        H2F hi_w_mm = new H2F("hi_w_mm", "hi_w_mm", 100, 0.6, ebeam*0.5+0.4, 100, 0.6, ebeam*0.4+0.6);
        hi_w_mm.setTitleX("W (GeV)");
        hi_w_mm.setTitleY("Mx2 (GeV)");
//        hi_dphi.setOptStat("1110");
        DataGroup dg_general = new DataGroup(2,3);
        dg_general.addDataSet(hi_rec_q2w,   0);
        dg_general.addDataSet(hi_rec_w,     1);
        dg_general.addDataSet(hi_rec_w_phi, 2);
        dg_general.addDataSet(hi_rec_el,    3);
        dg_general.addDataSet(hi_rec_mm,    4);  
        dg_general.addDataSet(f1_mm,        4);
        dg_general.addDataSet(hi_w_mm,      5);  
        this.getDataGroup().add(dg_general, 1);
        // MissingMass
        DataGroup dg_proton = new DataGroup(2,3);
        for(int sector=1; sector <= 6; sector++) {
            H1F hi_rec_mm_sec = new H1F("hi_rec_mm_" + sector, "hi_rec_mm_" + sector, 200, 0.6, ebeam*0.4+0.6);
            hi_rec_mm_sec.setTitleX("Mx (GeV)");
            hi_rec_mm_sec.setTitleY("Counts");
            hi_rec_mm_sec.setTitle("Sector " + sector); 
            F1D f1_mm_sec = new F1D("f1_mm_" + sector, "[amp]*gaus(x,[mean],[sigma])", 0.8, 1.2);
            f1_mm_sec.setParameter(0, 0);
            f1_mm_sec.setParameter(1, 1);
            f1_mm_sec.setParameter(2, 0.2);
            f1_mm_sec.setLineWidth(2);
            f1_mm_sec.setLineColor(2);
            f1_mm_sec.setOptStat("1111");
            dg_proton.addDataSet(hi_rec_mm_sec, sector-1);  
            dg_proton.addDataSet(f1_mm_sec    , sector-1);
        }
        this.getDataGroup().add(dg_proton, 2);   
    }
        
    @Override
    public void plotHistos() {
        this.getAnalysisCanvas().getCanvas("General").divide(3,2);
        this.getAnalysisCanvas().getCanvas("General").setGridX(false);
        this.getAnalysisCanvas().getCanvas("General").setGridY(false);
        this.getAnalysisCanvas().getCanvas("MissingMass").divide(3, 2);
        this.getAnalysisCanvas().getCanvas("MissingMass").setGridX(false);
        this.getAnalysisCanvas().getCanvas("MissingMass").setGridY(false);
        this.getAnalysisCanvas().getCanvas("General").cd(0);
        this.getAnalysisCanvas().getCanvas("General").draw(this.getDataGroup().getItem(1).getH2F("hi_rec_q2w"));
        this.getAnalysisCanvas().getCanvas("General").getPad(0).getAxisZ().setLog(true);        
        this.getAnalysisCanvas().getCanvas("General").cd(1);
        this.getAnalysisCanvas().getCanvas("General").draw(this.getDataGroup().getItem(1).getH1F("hi_rec_w"));
        this.getAnalysisCanvas().getCanvas("General").cd(2);
        this.getAnalysisCanvas().getCanvas("General").draw(this.getDataGroup().getItem(1).getH2F("hi_rec_w_phi"));
        this.getAnalysisCanvas().getCanvas("General").getPad(2).getAxisZ().setLog(true);        
        this.getAnalysisCanvas().getCanvas("General").cd(3);
        this.getAnalysisCanvas().getCanvas("General").draw(this.getDataGroup().getItem(1).getH2F("hi_rec_el"));
        this.getAnalysisCanvas().getCanvas("General").getPad(3).getAxisZ().setLog(true);        
        this.getAnalysisCanvas().getCanvas("General").cd(4);
        this.getAnalysisCanvas().getCanvas("General").draw(this.getDataGroup().getItem(1).getH1F("hi_rec_mm"));
        this.getAnalysisCanvas().getCanvas("General").draw(this.getDataGroup().getItem(1).getF1D("f1_mm"),"same");       
        this.getAnalysisCanvas().getCanvas("General").cd(5);
        this.getAnalysisCanvas().getCanvas("General").draw(this.getDataGroup().getItem(1).getH2F("hi_w_mm"));
        this.getAnalysisCanvas().getCanvas("General").getPad(5).getAxisZ().setLog(true);        
        this.getAnalysisCanvas().getCanvas("General").update();
         for(int sector=1; sector <= 6; sector++) {
            this.getAnalysisCanvas().getCanvas("MissingMass").cd(sector-1);
            this.getAnalysisCanvas().getCanvas("MissingMass").draw(this.getDataGroup().getItem(2).getH1F("hi_rec_mm_" + sector));
            this.getAnalysisCanvas().getCanvas("MissingMass").draw(this.getDataGroup().getItem(2).getF1D("f1_mm_" + sector),"same");
        }
        this.getAnalysisCanvas().getCanvas("MissingMass").update();        
    }
    
    @Override
    public void processEvent(DataEvent event) {
        int run = 0;
        if(event.hasBank("RUN::config")) {
            run = event.getBank("RUN::config").getInt("run", 0);       
        }
        else {
            return;
        }
        if(run == 0) return;
        double ebeamRCDB = 10.6;
        if(run>4000) ebeamRCDB = (double) this.getCcdb().getRcdbConstant(run, "beam_energy").getValue()/1000.;
        if(ebeamRCDB == 0) {
            ebeamRCDB = 10.6;
        }
        if(ebeamRCDB!=ebeam) {
            ebeam=ebeamRCDB;
            this.resetEventListener();
        }
         // process event info and save into data group
        Particle recEl = null;
        Particle recPi = null;
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
                if(bank.getInt("pid", loop)==11 && recEl==null && (short) Math.abs(bank.getShort("status", loop))>=2000) {
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
//                if(bank.getInt("charge", loop)==1 && recPi==null && bank.getShort("status", loop)<4000) {
                if(bank.getInt("pid", loop)==211 && recPi==null && (short) Math.abs(bank.getShort("status", loop))<4000) {
                    recPi = new Particle(
                                          211,
                                          bank.getFloat("px", loop),
                                          bank.getFloat("py", loop),
                                          bank.getFloat("pz", loop),
                                          bank.getFloat("vx", loop),
                                          bank.getFloat("vy", loop),
                                          bank.getFloat("vz", loop));
                }
            }
            if(recEl != null && recPi != null && nCharged==2) {
//            System.out.println("Analyzed ");
                virtualPhoton = new LorentzVector(0., 0., ebeam, ebeam);
                virtualPhoton.sub(recEl.vector());
                hadronSystem = new LorentzVector(0., 0., ebeam, 0.9383+ebeam);
                hadronSystem.sub(recEl.vector());
                missingParticle = new LorentzVector(0., 0., ebeam, 0.9383+ebeam);
                missingParticle.sub(recEl.vector());
                missingParticle.sub(recPi.vector());
                int secEl = (int) recEl.getProperty("sector");
                if(Math.toDegrees(recEl.theta())>0){
                    this.getDataGroup().getItem(1).getH2F("hi_rec_q2w").fill(hadronSystem.mass(),-virtualPhoton.mass2());
                    this.getDataGroup().getItem(1).getH1F("hi_rec_w").fill(hadronSystem.mass());
                    this.getDataGroup().getItem(1).getH2F("hi_rec_w_phi").fill(Math.toDegrees(recEl.phi()), hadronSystem.mass());
                    this.getDataGroup().getItem(1).getH2F("hi_rec_el").fill(recEl.p(),Math.toDegrees(recEl.theta()));
                }
                this.getDataGroup().getItem(1).getH1F("hi_rec_mm").fill(missingParticle.mass());
//                if(hadronSystem.mass()>1.1 && hadronSystem.mass()<1.4) {
                    this.getDataGroup().getItem(1).getH2F("hi_w_mm").fill(hadronSystem.mass(),missingParticle.mass());
                    this.getDataGroup().getItem(2).getH1F("hi_rec_mm_" + secEl).fill(missingParticle.mass());
//                }
            }

        }
    }

    @Override
    public void timerUpdate() {
        this.fitW(this.getDataGroup().getItem(1).getH1F("hi_rec_mm"), this.getDataGroup().getItem(1).getF1D("f1_mm"));
//        for(int sector=1; sector <= 6; sector++) {
//            this.fitW(this.getDataGroup().getItem(2).getH1F("hi_rec_mm_" + sector), this.getDataGroup().getItem(2).getF1D("f1_mm_" + sector));
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
