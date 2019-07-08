/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.analysis;

import java.util.ArrayList;
import org.clas.viewer.AnalysisMonitor;
import org.jlab.clas.pdg.PhysicsConstants;
import org.jlab.clas.physics.LorentzVector;
import org.jlab.clas.physics.Particle;
import org.jlab.detector.base.DetectorType;
import org.jlab.detector.calib.utils.ConstantsManager;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.F1D;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.utils.groups.IndexedTable;

/**
 *
 * @author devita
 */
public class KINEmonitor extends AnalysisMonitor {
    
    private double ebeam = 10.6;
    
    public KINEmonitor(String name, ConstantsManager ccdb) {
        super(name, ccdb);
        this.setAnalysisTabNames("Performance", "General");
        this.init(false);
    }

    
    @Override
    public void createHistos() {
        // initialize canvas and create histograms
        this.setNumberOfEvents(0);
        double maxW = 1.8;
        if(ebeam>6) maxW=3.5;
        
        H1F summary = new H1F("summary","summary",6,1,7);
        summary.setTitleX("sector");
        summary.setTitleY("DC hits");
        summary.setFillColor(33);
        DataGroup sum = new DataGroup(1,1);
        sum.addDataSet(summary, 0);
        this.setAnalysisSummary(sum);
        // General
        H2F hi_rec_q2x = new H2F("hi_rec_q2x","hi_rec_q2x",400, 0., 1.2, 400, 0., 14.); 
        hi_rec_q2x.setTitleX("x");
        hi_rec_q2x.setTitleY("Q2 (GeV2)");
        H2F hi_rec_q2w = new H2F("hi_rec_q2w","hi_rec_q2w",400, 0., 4.5, 400, 0., 14.); 
        hi_rec_q2w.setTitleX("W (GeV)");
        hi_rec_q2w.setTitleY("Q2 (GeV2)");
        DataGroup dg_general = new DataGroup(2,1);
        dg_general.addDataSet(hi_rec_q2w,   0);
        dg_general.addDataSet(hi_rec_q2x,   1); 
        this.getDataGroup().add(dg_general, 1);
        // beta
        DataGroup dc_part = new DataGroup(2,2);
        H2F hi_beta_pos_ftof = new H2F("hi_beta_pos_ftof", "hi_beta_pos_ftof", 500, 0., 5., 500, 0., 1.5);  
        hi_beta_pos_ftof.setTitleX("p (GeV)"); 
        hi_beta_pos_ftof.setTitleY("#beta");
        H2F hi_beta_pos_ctof = new H2F("hi_beta_pos_ctof", "hi_beta_pos_ctof", 500, 0., 3., 500, 0.2, 1.5);  
        hi_beta_pos_ctof.setTitleX("p (GeV)"); 
        hi_beta_pos_ctof.setTitleY("#beta");
        F1D fpion = new F1D("fpion","x/sqrt(x*x+0.1396*0.1396)", 0.2, 5.0);
        F1D fkaon = new F1D("fkaon","x/sqrt(x*x+0.4935*0.4935)", 0.2, 5.0);
        F1D fprot = new F1D("fprot","x/sqrt(x*x+0.9383*0.9383)", 0.2, 5.0);
        F1D cpion = new F1D("cpion","x/sqrt(x*x+0.1396*0.1396)", 0.2, 3.0);
        F1D ckaon = new F1D("ckaon","x/sqrt(x*x+0.4935*0.4935)", 0.2, 3.0);
        F1D cprot = new F1D("cprot","x/sqrt(x*x+0.9383*0.9383)", 0.2, 3.0);
        H1F hi_pi0_mass_fd = new H1F("hi_pi0_mass_fd", "hi_pi0_mass_fd", 200, 0.,400.);  
        hi_pi0_mass_fd.setTitleX("p(GeV)"); 
        hi_pi0_mass_fd.setTitleY("Counts");
        F1D fpi0_fd = new F1D("fpi0_fd", "[amp]*gaus(x,[mean],[sigma])+[p0]+[p1]*x", 50.,200.);
        fpi0_fd.setParameter(0, 0.0);
        fpi0_fd.setParameter(1, 140.0);
        fpi0_fd.setParameter(2, 2.0);
        fpi0_fd.setParameter(3, 0.0);
        fpi0_fd.setParameter(4, 0.0);
        fpi0_fd.setLineWidth(2);
        fpi0_fd.setLineColor(2);
        fpi0_fd.setOptStat("1111111");        
        H1F hi_pi0_mass_ft = new H1F("hi_pi0_mass_ft", "hi_pi0_mass_ft", 200, 0.,400.);  
        hi_pi0_mass_ft.setTitleX("p(GeV)"); 
        hi_pi0_mass_ft.setTitleY("Counts");
        F1D fpi0_ft = new F1D("fpi0_ft", "[amp]*gaus(x,[mean],[sigma])+[p0]+[p1]*x", 80.,200.);
        fpi0_ft.setParameter(0, 0.0);
        fpi0_ft.setParameter(1, 140.0);
        fpi0_ft.setParameter(2, 2.0);
        fpi0_ft.setParameter(3, 0.0);
        fpi0_ft.setParameter(4, 0.0);
        fpi0_ft.setLineWidth(2);
        fpi0_fd.setLineColor(2);
        fpi0_ft.setOptStat("1111111");
        dc_part.addDataSet(hi_beta_pos_ftof, 0);
        dc_part.addDataSet(fpion,  0);
        dc_part.addDataSet(fkaon,  0);
        dc_part.addDataSet(fprot,  0);
        dc_part.addDataSet(hi_beta_pos_ctof, 1);
        dc_part.addDataSet(cpion,  1);
        dc_part.addDataSet(ckaon,  1);
        dc_part.addDataSet(cprot,  1);
        dc_part.addDataSet(hi_pi0_mass_fd, 2);
        dc_part.addDataSet(fpi0_fd,        2);
        dc_part.addDataSet(hi_pi0_mass_ft, 3);
        dc_part.addDataSet(fpi0_ft,        3);
        this.getDataGroup().add(dc_part, 2);  

    }
        
    @Override
    public void plotHistos() {
        this.getAnalysisCanvas().getCanvas("General").divide(2,1);
        this.getAnalysisCanvas().getCanvas("General").setGridX(false);
        this.getAnalysisCanvas().getCanvas("General").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Performance").divide(2,2);
        this.getAnalysisCanvas().getCanvas("Performance").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Performance").setGridY(false);
        this.getAnalysisCanvas().getCanvas("General").cd(0);
        this.getAnalysisCanvas().getCanvas("General").draw(this.getDataGroup().getItem(1).getH2F("hi_rec_q2x"));
        this.getAnalysisCanvas().getCanvas("General").getPad(0).getAxisZ().setLog(true);        
        this.getAnalysisCanvas().getCanvas("General").cd(1);
        this.getAnalysisCanvas().getCanvas("General").draw(this.getDataGroup().getItem(1).getH2F("hi_rec_q2w"));
        this.getAnalysisCanvas().getCanvas("General").getPad(1).getAxisZ().setLog(true);        
        this.getAnalysisCanvas().getCanvas("Performance").cd(0);
        this.getAnalysisCanvas().getCanvas("Performance").getPad(0).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("Performance").draw(this.getDataGroup().getItem(2).getH2F("hi_beta_pos_ftof"));
        this.getAnalysisCanvas().getCanvas("Performance").draw(this.getDataGroup().getItem(2).getF1D("fpion"),"same");
        this.getAnalysisCanvas().getCanvas("Performance").draw(this.getDataGroup().getItem(2).getF1D("fkaon"),"same");
        this.getAnalysisCanvas().getCanvas("Performance").draw(this.getDataGroup().getItem(2).getF1D("fprot"),"same");
        this.getAnalysisCanvas().getCanvas("Performance").cd(1);
        this.getAnalysisCanvas().getCanvas("Performance").getPad(1).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("Performance").draw(this.getDataGroup().getItem(2).getH2F("hi_beta_pos_ctof"));
        this.getAnalysisCanvas().getCanvas("Performance").draw(this.getDataGroup().getItem(2).getF1D("cpion"),"same");
        this.getAnalysisCanvas().getCanvas("Performance").draw(this.getDataGroup().getItem(2).getF1D("ckaon"),"same");
        this.getAnalysisCanvas().getCanvas("Performance").draw(this.getDataGroup().getItem(2).getF1D("cprot"),"same");
        this.getAnalysisCanvas().getCanvas("Performance").cd(2);
        this.getAnalysisCanvas().getCanvas("Performance").draw(this.getDataGroup().getItem(2).getH1F("hi_pi0_mass_fd"));
        this.getAnalysisCanvas().getCanvas("Performance").draw(this.getDataGroup().getItem(2).getF1D("fpi0_fd"),"same");
        this.getAnalysisCanvas().getCanvas("Performance").cd(3);
        this.getAnalysisCanvas().getCanvas("Performance").draw(this.getDataGroup().getItem(2).getH1F("hi_pi0_mass_ft"));
        this.getAnalysisCanvas().getCanvas("Performance").draw(this.getDataGroup().getItem(2).getF1D("fpi0_ft"),"same");
    }
    
    @Override
    public void processEvent(DataEvent event) {
        // process event info and save into data group        
        // event builder
        DataBank recRun    = null;
        DataBank recBankEB = null;
        DataBank recDeteEB = null;
        DataBank recCaloEB = null;
        DataBank recFTagEB = null;
        DataBank recTracEB = null;
        DataBank recEvenEB = null;
        if(event.hasBank("RUN::config"))            recRun      = event.getBank("RUN::config");
        if(event.hasBank("REC::Particle"))          recBankEB   = event.getBank("REC::Particle");
        if(event.hasBank("REC::Scintillator"))      recDeteEB   = event.getBank("REC::Scintillator");
        if(event.hasBank("REC::Calorimeter"))       recCaloEB   = event.getBank("REC::Calorimeter");
        if(event.hasBank("REC::ForwardTagger"))     recFTagEB   = event.getBank("REC::ForwardTagger");
        if(event.hasBank("REC::Track"))             recTracEB   = event.getBank("REC::Track");
        if(event.hasBank("REC::Event"))             recEvenEB   = event.getBank("REC::Event");
        int ev  = 0;
        int run = 0;
        if(recRun!=null) {
            run = event.getBank("RUN::config").getInt("run", 0);  
//            if(run==5990) run=5900;
        }
        else {
            return;
        }
        IndexedTable rfConfig = this.getCcdb().getConstants(run, "/calibration/eb/rf/config");
        double rfPeriod = rfConfig.getDoubleValue("clock", 1,1,1);
        double ebeamRCDB = (double) this.getCcdb().getRcdbConstant(run, "beam_energy").getValue()/1000.;
        if(ebeamRCDB == 0) {
            ebeamRCDB = 10.6;
        }
        if(ebeamRCDB!=ebeam) {
            ebeam=ebeamRCDB;
            this.resetEventListener();
        }
        Particle recEl = null;
        ArrayList<Particle> gammasFD = new ArrayList();
        ArrayList<Particle> gammasFT = new ArrayList();
        LorentzVector virtualPhoton  = null;
        LorentzVector hadronSystem   = null;
        LorentzVector virtualPhotonP = null;
        LorentzVector hadronSystemP  = null;
        if(recBankEB!=null && recEvenEB!=null ){
            gammasFD.clear();
            gammasFT.clear();
            double startTime = recEvenEB.getFloat("startTime", 0);
            double    rfTime = recEvenEB.getFloat("RFTime", 0);
            int rows = recBankEB.rows();
            for(int loop = 0; loop < rows; loop++){
                int pid    = recBankEB.getInt("pid", loop);
                int charge = recBankEB.getByte("charge", loop);
                float px   = recBankEB.getFloat("px", loop);
                float py   = recBankEB.getFloat("py", loop);
                float pz   = recBankEB.getFloat("pz", loop);
                float vx   = recBankEB.getFloat("vx", loop);
                float vy   = recBankEB.getFloat("vy", loop);
                float vz   = recBankEB.getFloat("vz", loop);
                float beta = recBankEB.getFloat("beta", loop);
                short status = (short) Math.abs(recBankEB.getShort("status", loop));
                if(pid==0) {
                    if(charge!=0) pid=211*charge;
                    else          pid=22;
                }
                Particle recParticle = new Particle(pid,px,py,pz,vx,vy,vz);
                recParticle.setProperty("beta", (double) beta);
                recParticle.setProperty("status", (double) status);
                if(charge!=0  && recDeteEB!=null ) {
                    for(int i=0; i<recDeteEB.rows(); i++) {
                        int pindex   = recDeteEB.getShort("pindex", i);
                        int sector   = recDeteEB.getByte("sector", i);
                        int layer    = recDeteEB.getByte("layer", i);
                        int detector = recDeteEB.getByte("detector", i);
                        double time  = recDeteEB.getFloat("time", i);
                        double path  = recDeteEB.getFloat("path", i);
                        if(pindex==loop && ((detector==DetectorType.FTOF.getDetectorId() && layer==2) || detector==DetectorType.CTOF.getDetectorId())){
                            recParticle.setProperty("time", (double) time);
                            recParticle.setProperty("path", (double) path);
                            double dt = (recParticle.getProperty("time") - recParticle.getProperty("path")/(PhysicsConstants.speedOfLight()) - rfTime);
                            dt = (dt +1000.5*rfPeriod)%rfPeriod-0.5*rfPeriod;
                            recParticle.setProperty("vertexTime", dt);
                        }
                    }
                }
                
                if(loop==0 && pid==11 && status>=2000) {
                    recEl = new Particle(recParticle); 
                    if(recParticle.hasProperty("time")) {
                        recEl.setProperty("beta", recParticle.getProperty("beta"));
                        recEl.setProperty("time", recParticle.getProperty("time"));
                        recEl.setProperty("path", recParticle.getProperty("path"));
                        recEl.setProperty("vertexTime", recParticle.getProperty("vertexTime"));
                        recEl.setProperty("status", recParticle.getProperty("status"));
                    }
                }
                if(charge>0 && recEl!=null && recParticle.hasProperty("time") && recEl.getProperty("vertexTime")<0.5 && Math.abs(recEl.vz()+3)<5) {
                    if(recParticle.getProperty("status")>4000) this.getDataGroup().getItem(2).getH2F("hi_beta_pos_ctof").fill(recParticle.p(),beta);
                    else                                       this.getDataGroup().getItem(2).getH2F("hi_beta_pos_ftof").fill(recParticle.p(),beta);

                }
                else if(charge==0 && recCaloEB!=null) {
                    double energy1=0;
                    double energy4=0;
                    for(int i=0; i<recCaloEB.rows(); i++) {
                        if(recCaloEB.getShort("pindex",i)==loop && recCaloEB.getByte("detector",i)==DetectorType.ECAL.getDetectorId()) {
                            if(recCaloEB.getByte("layer",i) == 1) energy1 = recCaloEB.getFloat("energy",i);
                            if(recCaloEB.getByte("layer",i) == 4) energy4 = recCaloEB.getFloat("energy",i);
                        }
                    }                    
                    if(status>2000 && status<3000 && (status-2000)<100 && (status-2000)>10 && energy1>0.05 && energy4>0) {
                        gammasFD.add(recParticle);
                    }
                    else if(status>=1000 && status<2000) {
                         gammasFT.add(recParticle);                    
                    }
                }
            }          
            // fill kinematics plots
            if(recEl != null) {
                virtualPhoton = new LorentzVector(0., 0., ebeam, ebeam);
                virtualPhoton.sub(recEl.vector());
                hadronSystem = new LorentzVector(0., 0., ebeam, 0.9383+ebeam);
                hadronSystem.sub(recEl.vector());
                double x=-virtualPhoton.mass2()/2/0.9383/virtualPhoton.e();
                if(Math.toDegrees(recEl.theta())>0 && recEl.p()>1){
                    this.getDataGroup().getItem(1).getH2F("hi_rec_q2w").fill(hadronSystem.mass(),-virtualPhoton.mass2());
                    this.getDataGroup().getItem(1).getH2F("hi_rec_q2x").fill(x,-virtualPhoton.mass2());
                }
            }
            // fill pi0 plots
            if (gammasFT.size() >= 2) {
                for (int i1 = 0; i1 < gammasFT.size(); i1++) {
                    for (int i2 = i1 + 1; i2 < gammasFT.size(); i2++) {
                        Particle partGamma1 = gammasFT.get(i1);
                        Particle partGamma2 = gammasFT.get(i2);
                        Particle partPi0 = new Particle();
                        partPi0.copy(partGamma1);
                        partPi0.combine(partGamma2, +1);
                        double invmass = Math.sqrt(partPi0.mass2());
                        double x = (partGamma1.p() - partGamma2.p()) / (partGamma1.p() + partGamma2.p());
                        double angle = Math.toDegrees(Math.acos(partGamma1.cosTheta(partGamma2)));
                        if (angle > 1.5) {
                            this.getDataGroup().getItem(2).getH1F("hi_pi0_mass_ft").fill(invmass * 1000);
                        }
                    }
                }
            }
            if (gammasFD.size() >= 2) {
                for (int i1 = 0; i1 < gammasFD.size(); i1++) {
                    for (int i2 = i1 + 1; i2 < gammasFD.size(); i2++) {
                        Particle partGamma1 = gammasFD.get(i1);
                        Particle partGamma2 = gammasFD.get(i2);
                        Particle partPi0 = new Particle();
                        partPi0.copy(partGamma1);
                        partPi0.combine(partGamma2, +1);
                        double invmass = Math.sqrt(partPi0.mass2());
                        double x = (partGamma1.p() - partGamma2.p()) / (partGamma1.p() + partGamma2.p());
                        double angle = Math.toDegrees(Math.acos(partGamma1.cosTheta(partGamma2)));
//                        if(Math.sqrt(partGamma1.p() + partGamma2.p()>0.1){
                        if(angle>3 && angle<17)  {
                            this.getDataGroup().getItem(2).getH1F("hi_pi0_mass_fd").fill(invmass * 1000);
                        }
                    }
                }
            }

        }
    }

    @Override
    public void timerUpdate() {
        this.fitpi0(this.getDataGroup().getItem(2).getH1F("hi_pi0_mass_fd"),this.getDataGroup().getItem(2).getF1D("fpi0_fd"));
        this.fitpi0(this.getDataGroup().getItem(2).getH1F("hi_pi0_mass_ft"),this.getDataGroup().getItem(2).getF1D("fpi0_ft"));
    }

    public void fitpi0(H1F h1p, F1D f1p) {
        // fit pi0 mass
        double hAmp  = h1p.getBinContent(h1p.getMaximumBin());
        double hMean = h1p.getAxis().getBinCenter(h1p.getMaximumBin());
        double hRMS  = 10; //ns
        f1p.setParameter(0, hAmp);
        f1p.setParLimits(0, hAmp*0.8, hAmp*1.2);
        f1p.setParameter(1, hMean);
        f1p.setParLimits(1, hMean-hRMS, hMean+hRMS);
        DataFitter.fit(f1p,h1p,"LQ");
        h1p.setFunction(null);   
    }
}
