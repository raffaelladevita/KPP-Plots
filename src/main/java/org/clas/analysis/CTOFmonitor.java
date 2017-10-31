/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.analysis;

import java.util.ArrayList;
import org.clas.viewer.AnalysisMonitor;
import org.jlab.clas.physics.Particle;
import org.jlab.geom.prim.Vector3D;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.fitter.ParallelSliceFitter;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.F1D;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;

/**
 *
 * @author devita
 */
public class CTOFmonitor extends AnalysisMonitor {
    

    public CTOFmonitor(String name) {
        super(name);
        this.setAnalysisTabNames("MIPs", "dE/dx", "RF","Beta and Mass");
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
        // dE/dx
        DataGroup dc_energy = new DataGroup(3,2);
        Integer nPaddle = 48;
        H1F hi_pos_en = new H1F("hi_pos_en" , "hi_pos_en", 100, 0., 25.);   
        hi_pos_en.setTitleX("dE/dx (MeV/cm)");
        hi_pos_en.setTitleY("Counts");
        H2F hi_pos_en_p = new H2F("hi_pos_en_p", "hi_pos_en_p", 100, 0, 2, 100, 0., 25.); 
        hi_pos_en_p.setTitleX("p (GeV)"); 
        hi_pos_en_p.setTitleY("dE/dx (MeV/cm)");
        H2F hi_pos_en_paddle = new H2F("hi_pos_en_paddle", "hi_pos_en_paddle", nPaddle, 1, nPaddle+1., 100, 0., 25.); 
        hi_pos_en_paddle.setTitleX("Counter"); 
        hi_pos_en_paddle.setTitleY("dE/dx (MeV/cm)");
        H1F hi_neg_en = new H1F("hi_neg_en", "hi_neg_en", 100, 0., 25.); 
        hi_neg_en.setTitleX("dE/dx (MeV/cm)");
        hi_neg_en.setTitleY("Counts");
        H2F hi_neg_en_p = new H2F("hi_neg_en_p", "hi_neg_en_p", 100, 0, 2, 100, 0., 25.);  
        hi_neg_en_p.setTitleX("p (GeV)"); 
        hi_neg_en_p.setTitleY("dE/dx (MeV/cm)");
        H2F hi_neg_en_paddle = new H2F("hi_neg_en_paddle", "hi_neg_en_paddle", nPaddle, 1, nPaddle+1., 100, 0., 25.); 
        hi_neg_en_paddle.setTitleX("Counter"); 
        hi_neg_en_paddle.setTitleY("dE/dx (MeV/cm)");
        dc_energy.addDataSet(hi_pos_en,       0);
        dc_energy.addDataSet(hi_pos_en_p,     1);
        dc_energy.addDataSet(hi_pos_en_paddle,2);
        dc_energy.addDataSet(hi_neg_en,       3);
        dc_energy.addDataSet(hi_neg_en_p,     4);
        dc_energy.addDataSet(hi_neg_en_paddle,5);
        this.getDataGroup().add(dc_energy, 1);
        // paddle info
        DataGroup dc_mips = new DataGroup(1,2);
        H2F hi_en_paddle = new H2F("hi_en_paddle", "hi_en_paddle", nPaddle, 1, nPaddle+1., 100, 1, 51);  
        hi_en_paddle.setTitleX("Counter"); 
        hi_en_paddle.setTitleY("Energy (MeV)");
        H2F hi_time_paddle = new H2F("hi_time_paddle", "hi_time_paddle", nPaddle, 1, nPaddle+1., 100, -20., 20.);  
        hi_time_paddle.setTitleX("Counter"); 
        hi_time_paddle.setTitleY("Up-Down Time (ns)");
        dc_mips.addDataSet(hi_en_paddle,  0);
        dc_mips.addDataSet(hi_time_paddle,1);        
        this.getDataGroup().add(dc_mips, 2);  
        // RF offsets
        DataGroup dc_rf = new DataGroup(2,2);
        H2F hi_rf_neg_paddle = new H2F("hi_rf_neg_paddle", "hi_rf_neg_paddle", nPaddle, 1, nPaddle+1., 100, -1., 1.);  
        hi_rf_neg_paddle.setTitleX("Counter");
        hi_rf_neg_paddle.setTitleY("Vertex Time (ns)"); 
        H2F hi_rf_pos_paddle = new H2F("hi_rf_pos_paddle", "hi_rf_pos_paddle", nPaddle, 1, nPaddle+1., 100, -1., 1.);  
        hi_rf_pos_paddle.setTitleX("Counter");
        hi_rf_pos_paddle.setTitleY("Vertex Time (ns)"); 
        GraphErrors g_rf_neg_paddle = new GraphErrors("g_rf_neg_paddle");
        g_rf_neg_paddle.setTitleX("RF offset (ns)"); 
        g_rf_neg_paddle.setTitleY("Counter");
        GraphErrors g_rf_pos_paddle = new GraphErrors("g_rf_pos_paddle");
        g_rf_pos_paddle.setTitleX("RF offset (ns)"); 
        g_rf_pos_paddle.setTitleY("Counter");
        for(int i=1; i<=nPaddle; i++) g_rf_pos_paddle.addPoint((double) i, 0, 0, 0);
        dc_rf.addDataSet(hi_rf_neg_paddle, 0);
        dc_rf.addDataSet(hi_rf_pos_paddle, 1);
        dc_rf.addDataSet(g_rf_neg_paddle,  2);
        dc_rf.addDataSet(g_rf_pos_paddle,  3);        
        this.getDataGroup().add(dc_rf, 3);   
        // beta
        DataGroup dc_beta = new DataGroup(2,2);
        H1F hi_mass_pos = new H1F("hi_mass_pos", "hi_mass_pos", 100, -1., 3.);  
        hi_mass_pos.setTitleX("Mass^2 (GeV)"); 
        hi_mass_pos.setTitleY("Counts");
        hi_mass_pos.setFillColor(32);
        H1F hi_mass_neg = new H1F("hi_mass_neg", "hi_mass_neg", 100, -1., 3.);  
        hi_mass_neg.setTitleX("Mass^2 (GeV)"); 
        hi_mass_neg.setTitleY("Counts");
        hi_mass_neg.setFillColor(34);
        H2F hi_beta_pos = new H2F("hi_beta_pos", "hi_beta_pos", 100, 0., 2., 100, 0., 1.4);  
        hi_beta_pos.setTitleX("p (GeV)"); 
        hi_beta_pos.setTitleY("#beta");
        H2F hi_beta_neg = new H2F("hi_beta_neg", "hi_beta_neg", 100, 0., 2., 100, 0., 1.4);  
        hi_beta_neg.setTitleX("p (GeV)"); 
        hi_beta_neg.setTitleY("#beta");
        dc_beta.addDataSet(hi_mass_pos,  0);
        dc_beta.addDataSet(hi_mass_neg,  1);
        dc_beta.addDataSet(hi_beta_pos,  2);
        dc_beta.addDataSet(hi_beta_neg,  3);
        this.getDataGroup().add(dc_beta, 4);  
    }
    
    @Override
    public void plotHistos() {
        this.getAnalysisCanvas().getCanvas("MIPs").divide(2,1);
        this.getAnalysisCanvas().getCanvas("MIPs").setGridX(false);
        this.getAnalysisCanvas().getCanvas("MIPs").setGridY(false);
        this.getAnalysisCanvas().getCanvas("RF").divide(2, 2);
        this.getAnalysisCanvas().getCanvas("RF").setGridX(false);
        this.getAnalysisCanvas().getCanvas("RF").setGridY(false);
        this.getAnalysisCanvas().getCanvas("dE/dx").divide(2, 2);
        this.getAnalysisCanvas().getCanvas("dE/dx").setGridX(false);
        this.getAnalysisCanvas().getCanvas("dE/dx").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Beta and Mass").divide(2, 2);
        this.getAnalysisCanvas().getCanvas("Beta and Mass").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Beta and Mass").setGridY(false);
         
        this.getAnalysisCanvas().getCanvas("MIPs").cd(0);
        this.getAnalysisCanvas().getCanvas("MIPs").getPad(0).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("MIPs").draw(this.getDataGroup().getItem(2).getH2F("hi_en_paddle"));
        this.getAnalysisCanvas().getCanvas("MIPs").cd(1);
        this.getAnalysisCanvas().getCanvas("MIPs").getPad(1).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("MIPs").draw(this.getDataGroup().getItem(2).getH2F("hi_time_paddle"));
        this.getAnalysisCanvas().getCanvas("RF").cd(0);
        this.getAnalysisCanvas().getCanvas("RF").draw(this.getDataGroup().getItem(3).getH2F("hi_rf_pos_paddle"));
        this.getAnalysisCanvas().getCanvas("RF").cd(1);
        this.getAnalysisCanvas().getCanvas("RF").draw(this.getDataGroup().getItem(3).getH2F("hi_rf_neg_paddle"));
        this.getAnalysisCanvas().getCanvas("RF").cd(2);
        this.getAnalysisCanvas().getCanvas("RF").getPad(2).getAxisY().setRange(0, 0.3);
        this.getAnalysisCanvas().getCanvas("RF").draw(this.getDataGroup().getItem(3).getGraph("g_rf_pos_paddle"));
        this.getAnalysisCanvas().getCanvas("RF").cd(3);
        this.getAnalysisCanvas().getCanvas("RF").getPad(3).getAxisY().setRange(0, 0.3);
        this.getAnalysisCanvas().getCanvas("RF").draw(this.getDataGroup().getItem(3).getGraph("g_rf_neg_paddle"));
        this.getAnalysisCanvas().getCanvas("dE/dx").cd(0);
        this.getAnalysisCanvas().getCanvas("dE/dx").draw(this.getDataGroup().getItem(1).getH1F("hi_pos_en"));
        this.getAnalysisCanvas().getCanvas("dE/dx").cd(1);
        this.getAnalysisCanvas().getCanvas("dE/dx").draw(this.getDataGroup().getItem(1).getH1F("hi_neg_en"));
        this.getAnalysisCanvas().getCanvas("dE/dx").cd(2);
        this.getAnalysisCanvas().getCanvas("dE/dx").getPad(2).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("dE/dx").draw(this.getDataGroup().getItem(1).getH2F("hi_pos_en_p"));
        this.getAnalysisCanvas().getCanvas("dE/dx").cd(3);
        this.getAnalysisCanvas().getCanvas("dE/dx").getPad(3).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("dE/dx").draw(this.getDataGroup().getItem(1).getH2F("hi_neg_en_p"));
        this.getAnalysisCanvas().getCanvas("Beta and Mass").cd(0);
        this.getAnalysisCanvas().getCanvas("Beta and Mass").getPad(0).getAxisY().setLog(true);
        this.getAnalysisCanvas().getCanvas("Beta and Mass").draw(this.getDataGroup().getItem(4).getH1F("hi_mass_pos"));
        this.getAnalysisCanvas().getCanvas("Beta and Mass").cd(1);
        this.getAnalysisCanvas().getCanvas("Beta and Mass").getPad(1).getAxisY().setLog(true);
        this.getAnalysisCanvas().getCanvas("Beta and Mass").draw(this.getDataGroup().getItem(4).getH1F("hi_mass_neg"));
        this.getAnalysisCanvas().getCanvas("Beta and Mass").cd(2);
        this.getAnalysisCanvas().getCanvas("Beta and Mass").getPad(2).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("Beta and Mass").draw(this.getDataGroup().getItem(4).getH2F("hi_beta_pos"));
        this.getAnalysisCanvas().getCanvas("Beta and Mass").cd(3);
        this.getAnalysisCanvas().getCanvas("Beta and Mass").getPad(3).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("Beta and Mass").draw(this.getDataGroup().getItem(4).getH2F("hi_beta_neg"));

        this.getAnalysisCanvas().getCanvas("MIPs").update();    
        this.getAnalysisCanvas().getCanvas("RF").update();    
        this.getAnalysisCanvas().getCanvas("dE/dx").update();    
        this.getAnalysisCanvas().getCanvas("Beta and Mass").update();    
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
        DataBank recCtofHits = null;
        DataBank recCtofRaws = null;
        DataBank recHBTTrack = null;
        DataBank ctofADC = null;
        DataBank ctofTDC = null;
        if(event.hasBank("RUN::config"))            recRun      = event.getBank("RUN::config");
        if(event.hasBank("REC::Particle"))          recBankEB   = event.getBank("REC::Particle");
        if(event.hasBank("REC::Scintillator"))      recDeteEB   = event.getBank("REC::Scintillator");
        if(event.hasBank("REC::Event"))             recEvenEB   = event.getBank("REC::Event");
        if(event.hasBank("RUN::rf"))                recRunRF    = event.getBank("RUN::rf");
        if(event.hasBank("CTOF::hits"))             recCtofHits = event.getBank("CTOF::hits");
        if(event.hasBank("CTOF::rawhits"))          recCtofRaws = event.getBank("CTOF::rawhits");
        if(event.hasBank("CTOF::adc"))              ctofADC     = event.getBank("CTOF::adc");
        if(event.hasBank("CTOF::tdc"))              ctofTDC     = event.getBank("CTOF::tdc");
        if(event.hasBank("CVTRec::Tracks"))         recHBTTrack = event.getBank("CVTRec::Tracks");
        int ev = recRun.getInt("event",0);
//        System.out.println(ev); 
//            if(ev==134 || ev==370) {System.out.println(ev); recBankEB.show(); recDeteEB.show();}

        if(recCtofHits!=null) {
            int nrows = recCtofHits.rows();
            for(int loop = 0; loop < nrows; loop++){
                int paddle    = recCtofHits.getShort("component", loop);
                int trk_id    = recCtofHits.getShort("trkID", loop)+1;
                double energy = recCtofHits.getFloat("energy", loop);
                double time   = recCtofHits.getFloat("time",loop);
           	double tx     = recCtofHits.getFloat("tx", loop);
                double ty     = recCtofHits.getFloat("ty", loop);
                double tz     = recCtofHits.getFloat("tz", loop);
                double path   = recCtofHits.getFloat("pathLength", loop);
                double dx     = recCtofHits.getFloat("pathLengthThruBar", loop);
//        	int adcId1    = recCtofHits.getShort("adc_idx1", loop);
//                int adcId2    = recCtofHits.getShort("adc_idx2", loop);
//                int tdcId1    = recCtofHits.getShort("tdc_idx1", loop);
//                int tdcId2    = recCtofHits.getShort("tdc_idx2", loop);
//                double adc1   = (double) ctofADC.getInt("ADC",adcId1);
//                double adc2   = (double) ctofADC.getInt("ADC",adcId2);
//                double adct1  = (double) ctofADC.getFloat("time",adcId1);
//                double adct2  = (double) ctofADC.getFloat("time",adcId2);
//                double tdc1   = (double) ctofTDC.getInt("TDC",tdcId1);
//                double tdc2   = (double) ctofTDC.getInt("TDC",tdcId2);
                this.getDataGroup().getItem(2).getH2F("hi_en_paddle").fill(paddle*1.,energy);
                if(trk_id!=-1 && energy>1.5 && recRunRF!=null) {
        	    int    q    = recHBTTrack.getInt("q",trk_id-1);
                    double p    = recHBTTrack.getFloat("p",trk_id-1);
                    double pt   = recHBTTrack.getFloat("pt",trk_id-1);
                    double phi0 = recHBTTrack.getFloat("phi0",trk_id-1);
                    Particle recParticle = new Particle(211,pt*Math.cos(phi0),pt*Math.sin(phi0),Math.sqrt(p*p-pt*pt),0,0,0);
                    double beta = recParticle.p()/Math.sqrt(recParticle.p()*recParticle.p()+recParticle.mass2());
                    double trf = 0;
                    for(int k = 0; k < recRunRF.rows(); k++){
                        if(recRunRF.getInt("id", k)==1) trf = recRunRF.getFloat("time",k);
                    }
                    if(recEvenEB!=null) {
                        double startTime = recEvenEB.getFloat("STTime", 0);
                        if(startTime>-100) {
                            double dt = time - path/(beta*29.97) - startTime;
                            double betaTof = path/(time-startTime)/29.97;
                            double mass2   = Math.pow(recParticle.p()/betaTof, 2)-recParticle.p()*recParticle.p();
                            if(q==1)  {
                                this.getDataGroup().getItem(1).getH1F("hi_pos_en").fill(energy/dx);
                                this.getDataGroup().getItem(1).getH2F("hi_pos_en_p").fill(recParticle.p(),energy/dx);
                                this.getDataGroup().getItem(1).getH2F("hi_pos_en_paddle").fill(paddle*1.,energy/dx);
                                if(recParticle.p()>0.1) this.getDataGroup().getItem(3).getH2F("hi_rf_pos_paddle").fill(paddle*1.,dt);
                                this.getDataGroup().getItem(4).getH2F("hi_beta_pos").fill(recParticle.p(),betaTof);
                                this.getDataGroup().getItem(4).getH1F("hi_mass_pos").fill(mass2);
                            }
                            if(q==-1) {
                                this.getDataGroup().getItem(1).getH1F("hi_neg_en").fill(energy/dx);
                                this.getDataGroup().getItem(1).getH2F("hi_neg_en_p").fill(recParticle.p(),energy/dx);
                                this.getDataGroup().getItem(1).getH2F("hi_neg_en_paddle").fill(paddle*1.,energy/dx);
                                if(recParticle.p()>0.1) this.getDataGroup().getItem(3).getH2F("hi_rf_neg_paddle").fill(paddle*1.,dt);
                                this.getDataGroup().getItem(4).getH2F("hi_beta_neg").fill(recParticle.p(),betaTof);
                                this.getDataGroup().getItem(4).getH1F("hi_mass_neg").fill(mass2);
                            }
                        }
                    }
                }
            }
        }
        if(recCtofRaws!=null) {
            int nrows = recCtofRaws.rows();
            for(int loop = 0; loop < nrows; loop++){
                int paddle   = recCtofRaws.getShort("component", loop);
                float tleft  = recCtofRaws.getFloat("time_up", loop);
                float tright = recCtofRaws.getFloat("time_down", loop);
                this.getDataGroup().getItem(2).getH2F("hi_time_paddle").fill(paddle*1.,tleft-tright);
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
        for(int layer=1; layer <= 3; layer++) {
            H2F         h2    = this.getDataGroup().getItem(3).getH2F("hi_rf_neg_paddle");
            GraphErrors meanX = this.getDataGroup().getItem(3).getGraph("g_rf_neg_paddle");
            meanX.reset();
            ArrayList<H1F> hslice = h2.getSlicesX();
            for(int i=0; i<hslice.size(); i++) {
                double  x = h2.getXAxis().getBinCenter(i);
                double ex = 0;
                double  y = hslice.get(i).getRMS();
                double ey = 0;
                double mean  = hslice.get(i).getDataX(hslice.get(i).getMaximumBin());
                double amp   = hslice.get(i).getBinContent(hslice.get(i).getMaximumBin());
                double sigma = hslice.get(i).getRMS();
                F1D f1 = new F1D("f1_dvz_neg","[amp]*gaus(x,[mean],[sigma])+[p0]", -1.0, 1.0);
                f1.setParameter(0, amp);
                f1.setParameter(1, mean);
                f1.setParameter(2, 0.1);
                f1.setParameter(3, 0);
                DataFitter.fit(f1, hslice.get(i), "Q"); //No options uses error for sigma 
                if(amp>50) meanX.addPoint(x, f1.getParameter(2), ex, f1.parameter(2).error());
            }
//            this.getAnalysisCanvas().getCanvas("RF negative").cd(layer+2);
//            this.getAnalysisCanvas().getCanvas("RF negative").draw(meanY);
//            this.getAnalysisCanvas().getCanvas("RF negative").update();
        }
        for(int layer=1; layer <= 3; layer++) {
            H2F         h2    = this.getDataGroup().getItem(3).getH2F("hi_rf_pos_paddle");
            GraphErrors meanX = this.getDataGroup().getItem(3).getGraph("g_rf_pos_paddle");
            meanX.reset();
            ArrayList<H1F> hslice = h2.getSlicesX();
            for(int i=0; i<hslice.size(); i++) {
                double  x = h2.getXAxis().getBinCenter(i);
                double ex = 0;
                double  y = hslice.get(i).getRMS();
                double ey = 0;
                double mean  = hslice.get(i).getDataX(hslice.get(i).getMaximumBin());
                double amp   = hslice.get(i).getBinContent(hslice.get(i).getMaximumBin());
                double sigma = hslice.get(i).getRMS();
                F1D f1 = new F1D("f1_dvz_neg","[amp]*gaus(x,[mean],[sigma])+[p0]", -1.0, 1.0);
                f1.setParameter(0, amp);
                f1.setParameter(1, mean);
                f1.setParameter(2, 0.1);
                f1.setParameter(3, 0);
                DataFitter.fit(f1, hslice.get(i), "Q"); //No options uses error for sigma 
                if(amp>50) meanX.addPoint(x, f1.getParameter(2), ex, f1.parameter(2).error());
            }
//            this.getAnalysisCanvas().getCanvas("RF positive").cd(layer+2);
//            this.getAnalysisCanvas().getCanvas("RF positive").draw(meanX);
//            this.getAnalysisCanvas().getCanvas("RF positive").update();
        
        }
    }
}
