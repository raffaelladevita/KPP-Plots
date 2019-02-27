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
public class ELASTICmonitor extends AnalysisMonitor {
    
    private double ebeam = 2.22193;
    
    public ELASTICmonitor(String name) {
        super(name);
        this.setAnalysisTabNames("Beam", "Phi", "Electron", "Proton", "W", "General");
        this.init(false);
    }

    
    @Override
    public void createHistos() {
        // initialize canvas and create histograms
        this.setNumberOfEvents(0);
        double maxW = ebeam*0.6+0.6;
        double maxQ2 = 1;
        if(ebeam>6) {
            maxW  = 3.5;
            maxQ2 = 3;
        }
        H1F summary = new H1F("summary","summary",6,1,7);
        summary.setTitleX("sector");
        summary.setTitleY("DC hits");
        summary.setFillColor(33);
        DataGroup sum = new DataGroup(1,1);
        sum.addDataSet(summary, 0);
        this.setAnalysisSummary(sum);
        // General
        H2F hi_rec_q2w = new H2F("hi_rec_q2w","hi_rec_q2w",100, 0.6, ebeam*0.6+0.6, 100, 0., maxQ2); 
        hi_rec_q2w.setTitleX("W (GeV)");
        hi_rec_q2w.setTitleY("Q2 (GeV2)");
        H1F hi_rec_w = new H1F("hi_rec_w","hi_rec_w",250, 0.6, maxW); 
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
        H2F hi_rec_w_phi = new H2F("hi_rec_w_phi","hi_rec_w_phi",100, -180., 180., 250, 0.6, maxW); 
        hi_rec_w_phi.setTitleX("#phi (deg)");
        hi_rec_w_phi.setTitleY("W (GeV)");
        H2F hi_rec_el = new H2F("hi_rec_el","hi_rec_el",100, 0.5, ebeam+0.5, 100, 0., 35.);
        hi_rec_el.setTitleX("p (GeV)");
        hi_rec_el.setTitleY("#theta (deg)");
        hi_rec_el.setTitle("Electron");
        F1D f1_el = new F1D("f1_el", "2*(180/3.14)*atan(sqrt(0.93832*([e0]-x)/2/[e0]/x))", ebeam*0.75, ebeam*0.99);
        f1_el.setParameter(0, ebeam);
        H2F hi_rec_pr = new H2F("hi_rec_pr","hi_rec_pr",100, 0.1, ebeam/2, 100, 30., 90.);
        hi_rec_pr.setTitleX("p (GeV)");
        hi_rec_pr.setTitleY("#theta (deg)");
        hi_rec_pr.setTitle("Proton");       
        F1D f1_pr = new F1D("f1_pr", "(180/3.14)*acos(([e0]*[e0]+x*x-pow(([e0]+0.93832-sqrt(x*x+0.9382*0.9382)),2))/2/[e0]/x)", ebeam*0.1, ebeam*0.5);
        f1_pr.setParameter(0, ebeam);
        H2F hi_phi = new H2F("hi_phi", "hi_phi", 200, -180, 180, 200, -180, 180);   
        hi_phi.setTitleX("El #phi (deg)");
        hi_phi.setTitleY("Pr #phi (deg)");
//        hi_dphi.setOptStat("1110");
        DataGroup dg_general = new DataGroup(2,3);
        dg_general.addDataSet(hi_rec_q2w,   0);
        dg_general.addDataSet(hi_rec_w,     1);
        dg_general.addDataSet(f1_w,         1);
        dg_general.addDataSet(hi_rec_w_phi, 2);
        dg_general.addDataSet(hi_rec_el,    3);
        dg_general.addDataSet(f1_el,        3);
        dg_general.addDataSet(hi_rec_pr,    4);  
        dg_general.addDataSet(f1_pr,        4);
        dg_general.addDataSet(hi_phi,       5);  
        this.getDataGroup().add(dg_general, 1);
        // W
        DataGroup dg_w = new DataGroup(2,3);
        for(int sector=1; sector <= 6; sector++) {
            H1F hi_rec_w_sec = new H1F("hi_rec_w_" + sector, "hi_rec_w_" + sector, 250, 0.6, maxW);  
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
            dg_w.addDataSet(hi_rec_w_sec, sector-1);
            dg_w.addDataSet(f1_w_sec    , sector-1);
        }
        this.getDataGroup().add(dg_w, 2);   
        // Proton
        DataGroup dg_proton = new DataGroup(2,3);
        for(int sector=1; sector <= 6; sector++) {
            H2F hi_rec_pr_sec = new H2F("hi_rec_pr_" + sector, "hi_rec_pr_" + sector, 100, 0.1, ebeam/2, 100, 30., 90.);
            hi_rec_pr_sec.setTitleX("p (GeV)");
            hi_rec_pr_sec.setTitleY("#theta (deg)");
            hi_rec_pr_sec.setTitle("Sector " + sector); 
            dg_proton.addDataSet(hi_rec_pr_sec, sector-1);  
        }
        this.getDataGroup().add(dg_proton, 3);  
        // Electron
        DataGroup dg_electron = new DataGroup(2,3);
        for(int sector=1; sector <= 6; sector++) {
            H2F hi_rec_el_sec = new H2F("hi_rec_el_" + sector, "hi_rec_el_" + sector, 100, ebeam*0.75, ebeam*0.99, 100, ebeam*0.75, ebeam*0.99);  
            hi_rec_el_sec.setTitleX("p (GeV)");
            hi_rec_el_sec.setTitleY("p_calc (GeV)");
            hi_rec_el_sec.setTitle("Sector " + sector);
            F1D f1_el_sec = new F1D("f1_el_" + sector, "x", ebeam*0.75, ebeam*0.99);
//            f1_el_sec.setParameter(0, 0);
            H2F hi_rec_de_theta_sec = new H2F("hi_rec_de_theta_" + sector, "hi_rec_de_theta_" + sector, 100, -0.2, 0.2, 100, 5, 20);  
            hi_rec_de_theta_sec.setTitleX("#Delta p (GeV)");
            hi_rec_de_theta_sec.setTitleY("#theta (deg)");
            hi_rec_de_theta_sec.setTitle("Sector " + sector);
            dg_electron.addDataSet(hi_rec_el_sec      , sector-1);
            dg_electron.addDataSet(f1_el_sec          , sector-1);
            dg_electron.addDataSet(hi_rec_de_theta_sec, sector-1);
        }
        this.getDataGroup().add(dg_electron, 4);   
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
        this.getDataGroup().add(dg_phi, 5);   
        // Beam
        DataGroup dg_beam = new DataGroup(2,3);
        for(int sector=1; sector <= 6; sector++) {
            H1F hi_beam_sec = new H1F("hi_beam_" + sector, "hi_beam_" + sector, 100, ebeam*0.75, ebeam*1.2);  
            hi_beam_sec.setTitleX("Beam Energy (GeV)");
            hi_beam_sec.setTitleY("Counts");
            hi_beam_sec.setTitle("Sector " + sector);
            F1D f1_beam_sec = new F1D("f1_beam_" + sector, "[amp]*gaus(x,[mean],[sigma])", ebeam*0.9, ebeam*1.1);
            f1_beam_sec.setParameter(0, 0);
            f1_beam_sec.setParameter(1, 2.1);
            f1_beam_sec.setParameter(2, 0.2);
            f1_beam_sec.setLineWidth(2);
            f1_beam_sec.setLineColor(2);
            f1_beam_sec.setOptStat("1111");
            dg_beam.addDataSet(hi_beam_sec, sector-1);
            dg_beam.addDataSet(f1_beam_sec, sector-1);
        }
        this.getDataGroup().add(dg_beam, 6);   
    }
        
    @Override
    public void plotHistos() {
        this.getAnalysisCanvas().getCanvas("General").divide(3,2);
        this.getAnalysisCanvas().getCanvas("General").setGridX(false);
        this.getAnalysisCanvas().getCanvas("General").setGridY(false);
        this.getAnalysisCanvas().getCanvas("W").divide(3, 2);
        this.getAnalysisCanvas().getCanvas("W").setGridX(false);
        this.getAnalysisCanvas().getCanvas("W").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Proton").divide(3, 2);
        this.getAnalysisCanvas().getCanvas("Proton").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Proton").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Electron").divide(3, 2);
        this.getAnalysisCanvas().getCanvas("Electron").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Electron").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Phi").divide(3, 2);
        this.getAnalysisCanvas().getCanvas("Phi").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Phi").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Beam").divide(3, 2);
        this.getAnalysisCanvas().getCanvas("Beam").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Beam").setGridY(false);
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
        this.getAnalysisCanvas().getCanvas("General").draw(this.getDataGroup().getItem(1).getH2F("hi_rec_pr"));
        this.getAnalysisCanvas().getCanvas("General").draw(this.getDataGroup().getItem(1).getF1D("f1_pr"),"same");
        this.getAnalysisCanvas().getCanvas("General").getPad(4).getAxisZ().setLog(true);        
        this.getAnalysisCanvas().getCanvas("General").cd(5);
        this.getAnalysisCanvas().getCanvas("General").draw(this.getDataGroup().getItem(1).getH2F("hi_phi"));
        this.getAnalysisCanvas().getCanvas("General").getPad(5).getAxisZ().setLog(true);        
        this.getAnalysisCanvas().getCanvas("General").update();
        for(int sector=1; sector <= 6; sector++) {
            this.getAnalysisCanvas().getCanvas("W").cd(sector-1);
            this.getAnalysisCanvas().getCanvas("W").draw(this.getDataGroup().getItem(2).getH1F("hi_rec_w_" + sector));
            this.getAnalysisCanvas().getCanvas("W").draw(this.getDataGroup().getItem(2).getF1D("f1_w_" + sector),"same");
        }
        this.getAnalysisCanvas().getCanvas("W").update();
         for(int sector=1; sector <= 6; sector++) {
            this.getAnalysisCanvas().getCanvas("Proton").cd(sector-1);
            this.getAnalysisCanvas().getCanvas("Proton").getPad(sector-1).getAxisZ().setLog(true);
            this.getAnalysisCanvas().getCanvas("Proton").draw(this.getDataGroup().getItem(3).getH2F("hi_rec_pr_" + sector));
            this.getAnalysisCanvas().getCanvas("Proton").draw(this.getDataGroup().getItem(1).getF1D("f1_pr"),"same");
        }
        for(int sector=1; sector <= 6; sector++) {
            this.getAnalysisCanvas().getCanvas("Electron").cd(sector-1);
            this.getAnalysisCanvas().getCanvas("Electron").getPad(sector-1).getAxisZ().setLog(true);
//            this.getAnalysisCanvas().getCanvas("Electron").draw(this.getDataGroup().getItem(4).getH2F("hi_rec_el_" + sector));
//            this.getAnalysisCanvas().getCanvas("Electron").draw(this.getDataGroup().getItem(4).getF1D("f1_el_" + sector),"same");
            this.getAnalysisCanvas().getCanvas("Electron").draw(this.getDataGroup().getItem(4).getH2F("hi_rec_de_theta_" + sector));
        }
        this.getAnalysisCanvas().getCanvas("Proton").update();        
        for(int sector=1; sector <= 6; sector++) {
            this.getAnalysisCanvas().getCanvas("Phi").cd(sector-1);
            this.getAnalysisCanvas().getCanvas("Phi").draw(this.getDataGroup().getItem(5).getH1F("hi_dphi_" + sector));
            this.getAnalysisCanvas().getCanvas("Phi").draw(this.getDataGroup().getItem(5).getF1D("f1_dphi_" + sector),"same");
        }
        this.getAnalysisCanvas().getCanvas("Phi").update();        
        for(int sector=1; sector <= 6; sector++) {
            this.getAnalysisCanvas().getCanvas("Beam").cd(sector-1);
            this.getAnalysisCanvas().getCanvas("Beam").draw(this.getDataGroup().getItem(6).getH1F("hi_beam_" + sector));
            this.getAnalysisCanvas().getCanvas("Beam").draw(this.getDataGroup().getItem(6).getF1D("f1_beam_" + sector),"same");
        }
        this.getAnalysisCanvas().getCanvas("Beam").update();        
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
        double ebeamOld=ebeam;
        if(run>2365 && run<=2597)      ebeam=2.217;
        else if(run>3028 && run<=3105) ebeam=6.424;
        else if(run>3105 && run<=3817) ebeam=10.594;
        else if(run>3817 && run<=3861) ebeam=6.424;
        else if(run>3861 && run<=5670) ebeam=10.594;
        else if(run>5671 && run<=5874) ebeam=7.546;
        else if(run>5874)              ebeam=6.535;
        if(ebeamOld!=ebeam) {
            ebeamOld=ebeam;
            this.resetEventListener();
        }
        // process event info and save into data group
        Particle recEl = null;
        Particle recPr = null;
        LorentzVector virtualPhoton  = null;
        LorentzVector hadronSystem   = null;
        LorentzVector virtualPhotonP = null;
        LorentzVector hadronSystemP  = null;
        if(event.hasBank("REC::Particle")==true && event.hasBank("REC::Track")){
            DataBank  bank  = event.getBank("REC::Particle");
            DataBank  track = event.getBank("REC::Track");
            int rows = bank.rows();
            for(int loop = 0; loop < rows; loop++){
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
                else if(bank.getInt("charge", loop)==1 && recPr==null && bank.getShort("status", loop)>=2000) {
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
            if(recEl != null) {
//            System.out.println("Analyzed ");
                virtualPhoton = new LorentzVector(0., 0., ebeam, ebeam);
                virtualPhoton.sub(recEl.vector());
                hadronSystem = new LorentzVector(0., 0., ebeam, 0.9383+ebeam);
                hadronSystem.sub(recEl.vector());
                int secEl = (int) recEl.getProperty("sector");
                if(Math.toDegrees(recEl.theta())>0){
                    this.getDataGroup().getItem(1).getH2F("hi_rec_q2w").fill(hadronSystem.mass(),-virtualPhoton.mass2());
                    this.getDataGroup().getItem(1).getH1F("hi_rec_w").fill(hadronSystem.mass());
                    this.getDataGroup().getItem(1).getH2F("hi_rec_w_phi").fill(Math.toDegrees(recEl.phi()), hadronSystem.mass());
                    this.getDataGroup().getItem(1).getH2F("hi_rec_el").fill(recEl.p(),Math.toDegrees(recEl.theta()));
                    this.getDataGroup().getItem(2).getH1F("hi_rec_w_" + secEl).fill(hadronSystem.mass());
                }
                if(recPr != null) {
//                    virtualPhotonP = new LorentzVector(0., 0., 0, -0.9383);
//                    virtualPhoton.add(recPr.vector());
//                    hadronSystem = new LorentzVector(0., 0., ebeam, 0.9383+ebeam);
//                    hadronSystem.sub(recEl.vector());
                    this.getDataGroup().getItem(1).getH2F("hi_rec_pr").fill(recPr.p(),Math.toDegrees(recPr.theta()));
                    this.getDataGroup().getItem(4).getH2F("hi_rec_el_" + secEl).fill(recEl.p(),ebeam+0.93832-recPr.e());

                    double phPr = Math.toDegrees(recPr.phi()); 
                    double phEl = Math.toDegrees(recEl.phi()); 
                    if(phPr < phEl) phPr +=360;

                    this.getDataGroup().getItem(1).getH2F("hi_phi").fill(Math.toDegrees(recEl.phi()),Math.toDegrees(recPr.phi()));
                    this.getDataGroup().getItem(5).getH1F("hi_dphi_" + secEl).fill(phPr-phEl);
                    if(recPr != null && Math.abs(phPr-phEl-180)<10) {
                        this.getDataGroup().getItem(4).getH2F("hi_rec_de_theta_" + secEl).fill(recEl.p()-(ebeam+0.93832-recPr.e()),Math.toDegrees(recEl.theta()));
                        this.getDataGroup().getItem(3).getH2F("hi_rec_pr_" + secEl).fill(recPr.p(),Math.toDegrees(recPr.theta()));
                        this.getDataGroup().getItem(6).getH1F("hi_beam_" + secEl).fill((-0.93832+recPr.e()+recEl.p()));
                    }
                }
            }

        }
    }

    @Override
    public void timerUpdate() {
        this.fitW(this.getDataGroup().getItem(1).getH1F("hi_rec_w"), this.getDataGroup().getItem(1).getF1D("f1_w"));
        for(int sector=1; sector <= 6; sector++) {
            this.fitW(this.getDataGroup().getItem(2).getH1F("hi_rec_w_" + sector), this.getDataGroup().getItem(2).getF1D("f1_w_" + sector));
            this.fitPhi(this.getDataGroup().getItem(5).getH1F("hi_dphi_" + sector), this.getDataGroup().getItem(5).getF1D("f1_dphi_" + sector));
            this.fitEbeam(this.getDataGroup().getItem(6).getH1F("hi_beam_" + sector), this.getDataGroup().getItem(6).getF1D("f1_beam_" + sector));
        }
    }
    
    public void fitW(H1F hiw,F1D f1w) {

        // get histogram maximum in the rane 0.8-1.2
        int i1=hiw.getXaxis().getBin(0.8);
        int i2=hiw.getXaxis().getBin(1.25);
        double hiMax=0;
        int    imax=i1;
        for(int i=i1; i<=i2; i++) {
            if(hiMax<hiw.getBinContent(i)) {
                imax=i;
                hiMax=hiw.getBinContent(i);
            }
        }           
        double mean = hiw.getDataX(imax); //hiw.getDataX(hiw.getMaximumBin());
        double amp  = hiMax;//hiw.getBinContent(hiw.getMaximumBin());
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
    public void fitEbeam(H1F hiw,F1D f1w) {

        double mean = hiw.getDataX(hiw.getMaximumBin());
        double amp = hiw.getBinContent(hiw.getMaximumBin());
        double sigma = 0.04;
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
        rmin = mean - 1.5 * Math.abs(sigma);
        f1w.setRange(rmin, rmax);
//        System.out.println(mean + " " + sigma + " " + rmin + " " + rmax);
        DataFitter.fit(f1w, hiw, "Q"); //No options uses error for sigma 
        hiw.setFunction(null);
    }
}
