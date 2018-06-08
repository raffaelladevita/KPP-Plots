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
import org.jlab.groot.math.F1D;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;

/**
 *
 * @author devita
 */
public class EBmonitor extends AnalysisMonitor {
    

    public EBmonitor(String name) {
        super(name);
        this.setAnalysisTabNames("Vertex time", "Beta");
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
   
        // Vertex time
        DataGroup dc_time = new DataGroup(2,2);
        H2F hi_time_pos_ftof = new H2F("hi_time_pos_ftof", "hi_time_pos_ftof", 100, 0., 5., 100, -30, 30);  
        hi_time_pos_ftof.setTitleX("p (GeV)"); 
        hi_time_pos_ftof.setTitleY("Vertex time (ns)");
        H2F hi_time_neg_ftof = new H2F("hi_time_neg_ftof", "hi_time_neg_ftof", 100, 0., 5., 100, -30, 30); 
        hi_time_neg_ftof.setTitleX("p (GeV)"); 
        hi_time_neg_ftof.setTitleY("Vertex time (ns)");
        H2F hi_time_pos_ctof = new H2F("hi_time_pos_ctof", "hi_time_pos_ctof", 100, 0., 3., 100, -30, 30); 
        hi_time_pos_ctof.setTitleX("p (GeV)"); 
        hi_time_pos_ctof.setTitleY("Vertex time (ns)");
        H2F hi_time_neg_ctof = new H2F("hi_time_neg_ctof", "hi_time_neg_ctof", 100, 0., 3., 100, -30, 30);  
        hi_time_neg_ctof.setTitleX("p (GeV)"); 
        hi_time_neg_ctof.setTitleY("Vertex time (ns)");
        dc_time.addDataSet(hi_time_pos_ftof,  0);
        dc_time.addDataSet(hi_time_neg_ftof,  1);
        dc_time.addDataSet(hi_time_pos_ctof,  2);
        dc_time.addDataSet(hi_time_neg_ctof,  3);
        this.getDataGroup().add(dc_time, 0);  
        // beta
        DataGroup dc_beta = new DataGroup(2,2);
        H2F hi_beta_pos_ftof = new H2F("hi_beta_pos_ftof", "hi_beta_pos_ftof", 100, 0., 5., 100, 0., 1.5);  
        hi_beta_pos_ftof.setTitleX("p (GeV)"); 
        hi_beta_pos_ftof.setTitleY("#beta");
        H2F hi_beta_neg_ftof = new H2F("hi_beta_neg_ftof", "hi_beta_neg_ftof", 100, 0., 5., 100, 0., 1.5);  
        hi_beta_neg_ftof.setTitleX("p (GeV)"); 
        hi_beta_neg_ftof.setTitleY("#beta");
        H2F hi_beta_pos_ctof = new H2F("hi_beta_pos_ctof", "hi_beta_pos_ctof", 100, 0., 3., 100, 0., 1.5);  
        hi_beta_pos_ctof.setTitleX("p (GeV)"); 
        hi_beta_pos_ctof.setTitleY("#beta");
        H2F hi_beta_neg_ctof = new H2F("hi_beta_neg_ctof", "hi_beta_neg_ctof", 100, 0., 3., 100, 0., 1.5);  
        hi_beta_neg_ctof.setTitleX("p (GeV)"); 
        hi_beta_neg_ctof.setTitleY("#beta");
        F1D fpion = new F1D("fpion","x/sqrt(x*x+0.1396*0.1396)", 0.2, 5.0);
        F1D fkaon = new F1D("fkaon","x/sqrt(x*x+0.4935*0.4935)", 0.2, 5.0);
        F1D fprot = new F1D("fprot","x/sqrt(x*x+0.9383*0.9383)", 0.2, 5.0);
        F1D cpion = new F1D("cpion","x/sqrt(x*x+0.1396*0.1396)", 0.2, 3.0);
        F1D ckaon = new F1D("ckaon","x/sqrt(x*x+0.4935*0.4935)", 0.2, 3.0);
        F1D cprot = new F1D("cprot","x/sqrt(x*x+0.9383*0.9383)", 0.2, 3.0);
        dc_beta.addDataSet(hi_beta_pos_ftof,  0);
        dc_beta.addDataSet(fpion,  0);
        dc_beta.addDataSet(fkaon,  0);
        dc_beta.addDataSet(fprot,  0);
        dc_beta.addDataSet(hi_beta_neg_ftof,  1);
        dc_beta.addDataSet(fpion,  1);
        dc_beta.addDataSet(fkaon,  1);
        dc_beta.addDataSet(fprot,  1);
        dc_beta.addDataSet(hi_beta_pos_ctof,  2);
        dc_beta.addDataSet(cpion,  2);
        dc_beta.addDataSet(ckaon,  2);
        dc_beta.addDataSet(cprot,  2);
        dc_beta.addDataSet(hi_beta_neg_ctof,  3);
        dc_beta.addDataSet(cpion,  3);
        dc_beta.addDataSet(ckaon,  3);
        dc_beta.addDataSet(cprot,  3);        
        this.getDataGroup().add(dc_beta, 1);  
    }
    
    @Override
    public void plotHistos() {

        this.getAnalysisCanvas().getCanvas("Beta").divide(2, 2);
        this.getAnalysisCanvas().getCanvas("Beta").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Beta").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Vertex time").divide(2, 2);
        this.getAnalysisCanvas().getCanvas("Vertex time").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Vertex time").setGridY(false);
        
        
        this.getAnalysisCanvas().getCanvas("Beta").cd(0);
        this.getAnalysisCanvas().getCanvas("Beta").getPad(0).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("Beta").draw(this.getDataGroup().getItem(1).getH2F("hi_beta_pos_ftof"));
        this.getAnalysisCanvas().getCanvas("Beta").draw(this.getDataGroup().getItem(1).getF1D("fpion"),"same");
        this.getAnalysisCanvas().getCanvas("Beta").draw(this.getDataGroup().getItem(1).getF1D("fkaon"),"same");
        this.getAnalysisCanvas().getCanvas("Beta").draw(this.getDataGroup().getItem(1).getF1D("fprot"),"same");
        this.getAnalysisCanvas().getCanvas("Beta").cd(1);
        this.getAnalysisCanvas().getCanvas("Beta").getPad(1).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("Beta").draw(this.getDataGroup().getItem(1).getH2F("hi_beta_neg_ftof"));
        this.getAnalysisCanvas().getCanvas("Beta").draw(this.getDataGroup().getItem(1).getF1D("fpion"),"same");
        this.getAnalysisCanvas().getCanvas("Beta").draw(this.getDataGroup().getItem(1).getF1D("fkaon"),"same");
        this.getAnalysisCanvas().getCanvas("Beta").draw(this.getDataGroup().getItem(1).getF1D("fprot"),"same");
        this.getAnalysisCanvas().getCanvas("Beta").cd(2);
        this.getAnalysisCanvas().getCanvas("Beta").getPad(2).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("Beta").draw(this.getDataGroup().getItem(1).getH2F("hi_beta_pos_ctof"));
        this.getAnalysisCanvas().getCanvas("Beta").draw(this.getDataGroup().getItem(1).getF1D("cpion"),"same");
        this.getAnalysisCanvas().getCanvas("Beta").draw(this.getDataGroup().getItem(1).getF1D("ckaon"),"same");
        this.getAnalysisCanvas().getCanvas("Beta").draw(this.getDataGroup().getItem(1).getF1D("cprot"),"same");
        this.getAnalysisCanvas().getCanvas("Beta").cd(3);
        this.getAnalysisCanvas().getCanvas("Beta").getPad(3).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("Beta").draw(this.getDataGroup().getItem(1).getH2F("hi_beta_neg_ctof"));    
        this.getAnalysisCanvas().getCanvas("Beta").draw(this.getDataGroup().getItem(1).getF1D("cpion"),"same");
        this.getAnalysisCanvas().getCanvas("Beta").draw(this.getDataGroup().getItem(1).getF1D("ckaon"),"same");
        this.getAnalysisCanvas().getCanvas("Beta").draw(this.getDataGroup().getItem(1).getF1D("cprot"),"same");
        this.getAnalysisCanvas().getCanvas("Vertex time").cd(0);
        this.getAnalysisCanvas().getCanvas("Vertex time").getPad(0).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("Vertex time").draw(this.getDataGroup().getItem(0).getH2F("hi_time_pos_ftof"));
        this.getAnalysisCanvas().getCanvas("Vertex time").cd(1);
        this.getAnalysisCanvas().getCanvas("Vertex time").getPad(1).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("Vertex time").draw(this.getDataGroup().getItem(0).getH2F("hi_time_neg_ftof"));
        this.getAnalysisCanvas().getCanvas("Vertex time").cd(2);
        this.getAnalysisCanvas().getCanvas("Vertex time").getPad(2).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("Vertex time").draw(this.getDataGroup().getItem(0).getH2F("hi_time_pos_ctof"));
        this.getAnalysisCanvas().getCanvas("Vertex time").cd(3);
        this.getAnalysisCanvas().getCanvas("Vertex time").getPad(3).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("Vertex time").draw(this.getDataGroup().getItem(0).getH2F("hi_time_neg_ctof"));
    
        this.getAnalysisCanvas().getCanvas("Beta").update();    
        this.getAnalysisCanvas().getCanvas("Vertex time").update();    
    }
        
    @Override
    public void processEvent(DataEvent event) {
        // process event info and save into data group        
        // event builder
        DataBank recRun    = null;
        DataBank recBankEB = null;
        DataBank recDeteEB = null;
        DataBank recEvenEB = null;
        DataBank recRunRF  = null;
        if(event.hasBank("RUN::config"))            recRun      = event.getBank("RUN::config");
        if(event.hasBank("REC::Particle"))          recBankEB   = event.getBank("REC::Particle");
        if(event.hasBank("REC::Scintillator"))      recDeteEB   = event.getBank("REC::Scintillator");
        if(event.hasBank("REC::Event"))             recEvenEB   = event.getBank("REC::Event");
        if(event.hasBank("RUN::rf"))                recRunRF    = event.getBank("RUN::rf");
        int ev = recRun.getInt("event",0);
//        System.out.println(ev); 

        if(recEvenEB!=null && recBankEB!=null && recDeteEB!=null) {
            double startTime = recEvenEB.getFloat("STTime", 0);
            int nrows = recBankEB.rows();
            int trigger = 0;
            for(int i=0; i<nrows; i++) {
                int pid = recBankEB.getInt("pid", i);
                float px = recBankEB.getFloat("px", i);
                float py = recBankEB.getFloat("py", i);
                float pz = recBankEB.getFloat("pz", i);
                float vx = recBankEB.getFloat("vx", i);
                float vy = recBankEB.getFloat("vy", i);
                float vz = recBankEB.getFloat("vz", i);
                float beta = recBankEB.getFloat("beta", i);
                short status = recBankEB.getShort("status", i);
                Particle recParticle = null;
                if(pid!=0) {
                    recParticle = new Particle(pid,px,py,pz,vx,vy,vz);
                    recParticle.setProperty("status", (double) status);
                }
                if(i==0) trigger = pid;
                if(startTime>-100 && trigger==11 && i>0) {
                    
                    if(recParticle.getProperty("status")>=4000) {
                        if(recParticle.charge()==1)  {
                            this.getDataGroup().getItem(1).getH2F("hi_beta_pos_ctof").fill(recParticle.p(),beta);                    
                        }
                        if(recParticle.charge()==-1) {
                            this.getDataGroup().getItem(1).getH2F("hi_beta_neg_ctof").fill(recParticle.p(),beta);
                        }
                    }
                    else if(recParticle.getProperty("status")>=2000) {
                        if(recParticle.charge()==1)  {
                            this.getDataGroup().getItem(1).getH2F("hi_beta_pos_ftof").fill(recParticle.p(),beta);                    
                        }
                        if(recParticle.charge()==-1) {
                            this.getDataGroup().getItem(1).getH2F("hi_beta_neg_ftof").fill(recParticle.p(),beta);
                        }
                    }
                }
            }
            nrows = recDeteEB.rows();
            for(int i=0; i<nrows; i++) {
                int pindex   = recDeteEB.getShort("pindex", i);
                int detector = recDeteEB.getByte("detector", i);
                float time   = recDeteEB.getFloat("time", i);
                float path   = recDeteEB.getFloat("path", i);
                int charge   = recBankEB.getByte("charge", pindex);
                float px     = recBankEB.getFloat("px", pindex);
                float py     = recBankEB.getFloat("py", pindex);
                float pz     = recBankEB.getFloat("pz", pindex);
                float vx     = recBankEB.getFloat("vx", pindex);
                float vy     = recBankEB.getFloat("vy", pindex);
                float vz     = recBankEB.getFloat("vz", pindex);
                float beta   = recBankEB.getFloat("beta", pindex);
                int status   = recBankEB.getShort("status", pindex);
                if(charge!=0) {
                    Particle recParticle = new Particle(211*charge,px,py,pz,vx,vy,vz);
                    recParticle.setProperty("status", (double) status);
                    double betap = recParticle.p()/Math.sqrt(recParticle.p()*recParticle.p()+recParticle.mass2());
                    double dt = (time - path/(betap*29.97) - startTime);                

                    if(startTime>-100 && trigger==11 && pindex>0) {                    
                        if(recParticle.getProperty("status")>=4000) {
                            if(recParticle.charge()==1)  {
                                this.getDataGroup().getItem(0).getH2F("hi_time_pos_ctof").fill(recParticle.p(),dt);
                            }
                            if(recParticle.charge()==-1) {
                                this.getDataGroup().getItem(0).getH2F("hi_time_neg_ctof").fill(recParticle.p(),dt);
                            }
                        }
                        else if(recParticle.getProperty("status")>=2000) {
                            if(recParticle.charge()==1)  {
                                this.getDataGroup().getItem(0).getH2F("hi_time_pos_ftof").fill(recParticle.p(),dt);                    
                            }
                            if(recParticle.charge()==-1) {
                                this.getDataGroup().getItem(0).getH2F("hi_time_neg_ftof").fill(recParticle.p(),dt);
                            }
                        }
                    }
                }

            }

        }

   }

    @Override
    public void timerUpdate() {
//        System.out.println("Updating FTOF");
        this.analyze();
        //fitting negative tracks vertex
//        ParallelSliceFitter fitter = new ParallelSliceFitter(this.getDataGroup().getItem(4).getH2F("hi_rf_paddle_2"));
//        fitter.fitSlicesX();
//        GraphErrors meanY = fitter.getMeanSlices();
//        this.getDataGroup().getItem(4).getGraph("hi_rf_offsets_2").copy(meanY);       
    }
    
    @Override
    public void analyze() {
//        ParallelSliceFitter fitter = new ParallelSliceFitter(this.getDataGroup().getItem(4).getH2F("hi_rf_paddle_2"));
//        fitter.fitSlicesY();
//        GraphErrors meanY = fitter.getMeanSlices();
    }
}
