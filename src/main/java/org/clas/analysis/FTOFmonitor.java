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
public class FTOFmonitor extends AnalysisMonitor {
    

    public FTOFmonitor(String name) {
        super(name);
        this.setAnalysisTabNames("MIPs", "dE/dx positive", "dE/dx negative","RF positive","RF negative", "Beta", "Mass");
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
        // positive particles 
        DataGroup dg_positive = new DataGroup(3,2);
        String[] panels   = {"1A", "1B", "2"};
        Integer[] nPaddle = {23, 62, 5};
        for(int layer=1; layer <= 3; layer++) {
            H1F hi_pos_en = new H1F("hi_pos_en_" + layer, "hi_pos_en_" + layer, 100, 0.,10.);   
            hi_pos_en.setTitleX("dE/dx (MeV/cm)");
            hi_pos_en.setTitleY("Counts");
            hi_pos_en.setTitle("Panel " + panels[layer-1]);
            H2F hi_pos_en_p = new H2F("hi_pos_en_p_" + layer, "hi_pos_en_p_" + layer, 100, 0, 5, 100, 0, 10.);  
            hi_pos_en_p.setTitleX("p (GeV)"); 
            hi_pos_en_p.setTitleY("dE/dx (MeV/cm)");
            hi_pos_en_p.setTitle("Panel " + panels[layer-1]);
            H2F hi_pos_en_paddle = new H2F("hi_pos_en_paddle_" + layer, "hi_pos_en_paddle_" + layer, nPaddle[layer-1], 1, nPaddle[layer-1]+1., 100, 0, 10);  
            hi_pos_en_paddle.setTitleX("Counter"); 
            hi_pos_en_paddle.setTitleY("dE/dx (MeV/cm)");
            hi_pos_en_paddle.setTitle("Panel " + panels[layer-1]);
            dg_positive.addDataSet(hi_pos_en,      -1+layer);
            dg_positive.addDataSet(hi_pos_en_p,     2+layer);
            dg_positive.addDataSet(hi_pos_en_paddle,5+layer);
        }
        this.getDataGroup().add(dg_positive, 1);
        // negative particles 
        DataGroup dc_negative = new DataGroup(3,2);
        for(int layer=1; layer <= 3; layer++) {
            H1F hi_neg_en = new H1F("hi_neg_en_" + layer, "hi_neg_en_" + layer, 100, 0.,10.);   
            hi_neg_en.setTitleX("dE/dx (MeV/cm)");
            hi_neg_en.setTitleY("Counts");
            hi_neg_en.setTitle("Panel " + panels[layer-1]);
            H2F hi_neg_en_p = new H2F("hi_neg_en_p_" + layer, "hi_neg_en_p_" + layer, 100, 0, 5, 100, 0, 10.);  
            hi_neg_en_p.setTitleX("p (GeV)"); 
            hi_neg_en_p.setTitleY("dE/dx (MeV/cm)");
            hi_neg_en_p.setTitle("Panel " + panels[layer-1]);
            H2F hi_neg_en_paddle = new H2F("hi_neg_en_paddle_" + layer, "hi_neg_en_paddle_" + layer, nPaddle[layer-1], 1, nPaddle[layer-1]+1., 100, 0, 10);  
            hi_neg_en_paddle.setTitleX("Counter"); 
            hi_neg_en_paddle.setTitleY("dE/dx (MeV/cm)");
            hi_neg_en_paddle.setTitle("Panel " + panels[layer-1]);            
            dc_negative.addDataSet(hi_neg_en,      -1+layer);
            dc_negative.addDataSet(hi_neg_en_p,     2+layer);
            dc_negative.addDataSet(hi_neg_en_paddle,5+layer);
        }
        this.getDataGroup().add(dc_negative, 2);  
        // paddle info
        DataGroup dc_mips = new DataGroup(3,2);
        for(int layer=1; layer <= 3; layer++) {
            H2F hi_en_paddle = new H2F("hi_en_paddle_" + layer, "hi_en_paddle_" + layer, nPaddle[layer-1], 1, nPaddle[layer-1]+1., 100, 0, 30);  
            hi_en_paddle.setTitleX("Counter"); 
            hi_en_paddle.setTitleY("Energy (MeV)");
            hi_en_paddle.setTitle("Panel " + panels[layer-1]);
            H2F hi_time_paddle = new H2F("hi_time_paddle_" + layer, "hi_time_paddle_" + layer, 100, -40., 40., nPaddle[layer-1], 1, nPaddle[layer-1]+1.);  
            hi_time_paddle.setTitleX("Left-right Time (ns)"); 
            hi_time_paddle.setTitleY("Counter");
            hi_time_paddle.setTitle("Panel " + panels[layer-1]);
            dc_mips.addDataSet(hi_en_paddle,  -1+layer);
            dc_mips.addDataSet(hi_time_paddle, 2+layer);
        }
        this.getDataGroup().add(dc_mips, 3);  
        // RF offsets
        DataGroup dc_rf = new DataGroup(3,2);
        for(int layer=1; layer <= 3; layer++) {
            H2F hi_rf_paddle = new H2F("hi_rf_neg_paddle_" + layer, "hi_rf_neg_paddle_" + layer, nPaddle[layer-1], 1, nPaddle[layer-1]+1., 100, -1., 1.);  
            hi_rf_paddle.setTitleX("Counter");
            hi_rf_paddle.setTitleY("RF offset (ns)"); 
            hi_rf_paddle.setTitle("Panel " + panels[layer-1]);
            GraphErrors g_rf_paddle = new GraphErrors("g_rf_neg_paddle_" + layer);
            g_rf_paddle.setTitleX("RF offset (ns)"); 
            g_rf_paddle.setTitleY("Counter");
            g_rf_paddle.setTitle("Panel " + panels[layer-1]);  
            for(int i=1; i<=nPaddle[layer-1]; i++) g_rf_paddle.addPoint((double) i, 0, 0, 0);
            dc_rf.addDataSet(hi_rf_paddle, -1+layer);
            dc_rf.addDataSet(g_rf_paddle,   2+layer);
        }
        this.getDataGroup().add(dc_rf, 4);  
        // RF offsets hadrons
        DataGroup dc_rf_had = new DataGroup(3,2);
        for(int layer=1; layer <= 3; layer++) {
            H2F hi_rf_paddle = new H2F("hi_rf_pos_paddle_" + layer, "hi_rf_pos_paddle_" + layer, nPaddle[layer-1], 1, nPaddle[layer-1]+1., 100, -1., 1.);  
            hi_rf_paddle.setTitleX("Counter");
            hi_rf_paddle.setTitleY("RF offset (ns)"); 
            hi_rf_paddle.setTitle("Panel " + panels[layer-1]);
            GraphErrors g_rf_paddle = new GraphErrors("g_rf_pos_paddle_" + layer);
            g_rf_paddle.setTitleX("RF offset (ns)"); 
            g_rf_paddle.setTitleY("Counter");
            g_rf_paddle.setTitle("Panel " + panels[layer-1]);  
            for(int i=1; i<=nPaddle[layer-1]; i++) g_rf_paddle.addPoint((double) i, 0, 0, 0);
            dc_rf_had.addDataSet(hi_rf_paddle, -1+layer);
            dc_rf_had.addDataSet(g_rf_paddle,   2+layer);
        }
        this.getDataGroup().add(dc_rf_had, 5);  
        // beta
        DataGroup dc_beta = new DataGroup(3,2);
        for(int layer=1; layer <= 3; layer++) {
            H2F hi_beta_pos = new H2F("hi_beta_pos_" + layer, "hi_beta_pos_" + layer, 100, 0., 6., 100, 0., 1.2);  
            hi_beta_pos.setTitleX("p (GeV)"); 
            hi_beta_pos.setTitleY("#beta");
            hi_beta_pos.setTitle("Panel " + panels[layer-1]);
            H2F hi_beta_neg = new H2F("hi_beta_neg_" + layer, "hi_beta_neg_" + layer, 100, 0., 6., 100, 0., 1.2);  
            hi_beta_neg.setTitleX("p (GeV)"); 
            hi_beta_neg.setTitleY("#beta");
            hi_beta_neg.setTitle("Panel " + panels[layer-1]);
            dc_beta.addDataSet(hi_beta_pos,  -1+layer);
            dc_beta.addDataSet(hi_beta_neg,   2+layer);
        }
        this.getDataGroup().add(dc_beta, 6);  
        // beta
        DataGroup dc_mass = new DataGroup(3,2);
        for(int layer=1; layer <= 3; layer++) {
            H1F hi_mass_pos = new H1F("hi_mass_pos_" + layer, "hi_mass_pos_" + layer, 100, -1., 3.);  
            hi_mass_pos.setTitleX("Mass^2 (GeV)"); 
            hi_mass_pos.setTitleY("Counts");
            hi_mass_pos.setTitle("Panel " + panels[layer-1]);
            hi_mass_pos.setFillColor(32);
            H1F hi_mass_neg = new H1F("hi_mass_neg_" + layer, "hi_mass_neg_" + layer, 100, -1., 3.);  
            hi_mass_neg.setTitleX("Mass^2 (GeV)"); 
            hi_mass_neg.setTitleY("Counts");
            hi_mass_neg.setTitle("Panel " + panels[layer-1]);
            hi_mass_neg.setFillColor(34);
            dc_beta.addDataSet(hi_mass_pos,  -1+layer);
            dc_beta.addDataSet(hi_mass_neg,   2+layer);
        }
        this.getDataGroup().add(dc_beta, 7);  
    }
    
    @Override
    public void plotHistos() {
        this.getAnalysisCanvas().getCanvas("MIPs").divide(3, 2);
        this.getAnalysisCanvas().getCanvas("MIPs").setGridX(false);
        this.getAnalysisCanvas().getCanvas("MIPs").setGridY(false);
        this.getAnalysisCanvas().getCanvas("RF negative").divide(3, 2);
        this.getAnalysisCanvas().getCanvas("RF negative").setGridX(false);
        this.getAnalysisCanvas().getCanvas("RF negative").setGridY(false);
        this.getAnalysisCanvas().getCanvas("RF positive").divide(3, 2);
        this.getAnalysisCanvas().getCanvas("RF positive").setGridX(false);
        this.getAnalysisCanvas().getCanvas("RF positive").setGridY(false);
        this.getAnalysisCanvas().getCanvas("dE/dx positive").divide(3, 2);
        this.getAnalysisCanvas().getCanvas("dE/dx positive").setGridX(false);
        this.getAnalysisCanvas().getCanvas("dE/dx positive").setGridY(false);
        this.getAnalysisCanvas().getCanvas("dE/dx negative").divide(3, 2);
        this.getAnalysisCanvas().getCanvas("dE/dx negative").setGridX(false);
        this.getAnalysisCanvas().getCanvas("dE/dx negative").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Beta").divide(3, 2);
        this.getAnalysisCanvas().getCanvas("Beta").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Beta").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Mass").divide(3, 2);
        this.getAnalysisCanvas().getCanvas("Mass").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Mass").setGridY(false);
         
        for(int layer=1; layer <= 3; layer++) {
            this.getAnalysisCanvas().getCanvas("MIPs").cd(0+layer-1);
            this.getAnalysisCanvas().getCanvas("MIPs").getPad(0+layer-1).getAxisZ().setLog(true);
            this.getAnalysisCanvas().getCanvas("MIPs").draw(this.getDataGroup().getItem(3).getH2F("hi_en_paddle_"+layer));
            this.getAnalysisCanvas().getCanvas("MIPs").cd(3+layer-1);
            this.getAnalysisCanvas().getCanvas("MIPs").getPad(3+layer-1).getAxisZ().setLog(true);
            this.getAnalysisCanvas().getCanvas("MIPs").draw(this.getDataGroup().getItem(3).getH2F("hi_time_paddle_"+layer));
        }
        for(int layer=1; layer <= 3; layer++) {
            this.getAnalysisCanvas().getCanvas("RF positive").cd(0+layer-1);
            this.getAnalysisCanvas().getCanvas("RF positive").draw(this.getDataGroup().getItem(5).getH2F("hi_rf_pos_paddle_"+layer));
            this.getAnalysisCanvas().getCanvas("RF positive").cd(3+layer-1);
            this.getAnalysisCanvas().getCanvas("RF positive").getPad(3+layer-1).getAxisY().setRange(0, 0.3);
            this.getAnalysisCanvas().getCanvas("RF positive").draw(this.getDataGroup().getItem(5).getGraph("g_rf_pos_paddle_"+layer));
        }
        for(int layer=1; layer <= 3; layer++) {
            this.getAnalysisCanvas().getCanvas("RF negative").cd(0+layer-1);
            this.getAnalysisCanvas().getCanvas("RF negative").getPad(0+layer-1).getAxisZ().setLog(true);
            this.getAnalysisCanvas().getCanvas("RF negative").draw(this.getDataGroup().getItem(4).getH2F("hi_rf_neg_paddle_"+layer));
            this.getAnalysisCanvas().getCanvas("RF negative").cd(3+layer-1);
            this.getAnalysisCanvas().getCanvas("RF negative").getPad(3+layer-1).getAxisY().setRange(0, 0.3);
            this.getAnalysisCanvas().getCanvas("RF negative").draw(this.getDataGroup().getItem(4).getGraph("g_rf_neg_paddle_"+layer));
        }
        for(int layer=1; layer <= 3; layer++) {
            this.getAnalysisCanvas().getCanvas("dE/dx positive").cd(0+layer-1);
            this.getAnalysisCanvas().getCanvas("dE/dx positive").draw(this.getDataGroup().getItem(1).getH1F("hi_pos_en_"+layer));
            this.getAnalysisCanvas().getCanvas("dE/dx positive").cd(3+layer-1);
            this.getAnalysisCanvas().getCanvas("dE/dx positive").getPad(3+layer-1).getAxisZ().setLog(true);
            this.getAnalysisCanvas().getCanvas("dE/dx positive").draw(this.getDataGroup().getItem(1).getH2F("hi_pos_en_p_"+layer));
        }
        for(int layer=1; layer <= 3; layer++) {
            this.getAnalysisCanvas().getCanvas("dE/dx negative").cd(0+layer-1);
            this.getAnalysisCanvas().getCanvas("dE/dx negative").draw(this.getDataGroup().getItem(2).getH1F("hi_neg_en_"+layer));
            this.getAnalysisCanvas().getCanvas("dE/dx negative").cd(3+layer-1);
            this.getAnalysisCanvas().getCanvas("dE/dx negative").getPad(3+layer-1).getAxisZ().setLog(true);
            this.getAnalysisCanvas().getCanvas("dE/dx negative").draw(this.getDataGroup().getItem(2).getH2F("hi_neg_en_p_"+layer));
        }
        for(int layer=1; layer <= 3; layer++) {
            this.getAnalysisCanvas().getCanvas("Beta").cd(0+layer-1);
            this.getAnalysisCanvas().getCanvas("Beta").getPad(0+layer-1).getAxisZ().setLog(true);
            this.getAnalysisCanvas().getCanvas("Beta").draw(this.getDataGroup().getItem(6).getH2F("hi_beta_pos_"+layer));
            this.getAnalysisCanvas().getCanvas("Beta").cd(3+layer-1);
            this.getAnalysisCanvas().getCanvas("Beta").getPad(3+layer-1).getAxisZ().setLog(true);
            this.getAnalysisCanvas().getCanvas("Beta").draw(this.getDataGroup().getItem(6).getH2F("hi_beta_neg_"+layer));
        }
        for(int layer=1; layer <= 3; layer++) {
            this.getAnalysisCanvas().getCanvas("Mass").cd(0+layer-1);
            this.getAnalysisCanvas().getCanvas("Mass").getPad(0+layer-1).getAxisY().setLog(true);
            this.getAnalysisCanvas().getCanvas("Mass").draw(this.getDataGroup().getItem(7).getH1F("hi_mass_pos_"+layer));
            this.getAnalysisCanvas().getCanvas("Mass").cd(3+layer-1);
            this.getAnalysisCanvas().getCanvas("Mass").getPad(3+layer-1).getAxisY().setLog(true);
            this.getAnalysisCanvas().getCanvas("Mass").draw(this.getDataGroup().getItem(7).getH1F("hi_mass_neg_"+layer));
        }

        this.getAnalysisCanvas().getCanvas("MIPs").update();    
        this.getAnalysisCanvas().getCanvas("RF positive").update();    
        this.getAnalysisCanvas().getCanvas("RF negative").update();    
        this.getAnalysisCanvas().getCanvas("dE/dx positive").update();    
        this.getAnalysisCanvas().getCanvas("dE/dx negative").update();    
        this.getAnalysisCanvas().getCanvas("Beta").update();    
        this.getAnalysisCanvas().getCanvas("Mass").update();    
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
        DataBank recFtofHits = null;
        DataBank recFtofRaws = null;
        DataBank recHBTTrack = null;
        DataBank ftofADC = null;
        DataBank ftofTDC = null;
        if(event.hasBank("RUN::config"))            recRun      = event.getBank("RUN::config");
        if(event.hasBank("REC::Particle"))          recBankEB   = event.getBank("REC::Particle");
        if(event.hasBank("REC::Scintillator"))      recDeteEB   = event.getBank("REC::Scintillator");
        if(event.hasBank("REC::Event"))             recEvenEB   = event.getBank("REC::Event");
        if(event.hasBank("RUN::rf"))                recRunRF    = event.getBank("RUN::rf");
        if(event.hasBank("FTOF::hits"))             recFtofHits = event.getBank("FTOF::hits");
        if(event.hasBank("FTOF::rawhits"))          recFtofRaws = event.getBank("FTOF::rawhits");
        if(event.hasBank("FTOF::adc"))              ftofADC     = event.getBank("FTOF::adc");
        if(event.hasBank("FTOF::tdc"))               ftofTDC     = event.getBank("FTOF::tdc");
        if(event.hasBank("TimeBasedTrkg::TBTracks")) recHBTTrack = event.getBank("TimeBasedTrkg::TBTracks");
        int ev = recRun.getInt("event",0);
//        System.out.println(ev); 
//            if(ev==134 || ev==370) {System.out.println(ev); recBankEB.show(); recDeteEB.show();}
//        if(recBankEB!=null && recDeteEB!=null) {
//            int nrows = recBankEB.rows();
//            for(int loop = 0; loop < nrows; loop++){
//                int pidCode = 0;
//                if(recBankEB.getInt("pid", loop)!=0) pidCode = recBankEB.getInt("pid", loop);
//                else if(recBankEB.getByte("charge", loop)==-1) pidCode = -211;
//                else if(recBankEB.getByte("charge", loop)==1) pidCode = 211;
//                else pidCode = 2112;
//                Particle recParticle = new Particle(
//                                            pidCode,
//                                            recBankEB.getFloat("px", loop),
//                                            recBankEB.getFloat("py", loop),
//                                            recBankEB.getFloat("pz", loop),
//                                            recBankEB.getFloat("vx", loop),
//                                            recBankEB.getFloat("vy", loop),
//                                            recBankEB.getFloat("vz", loop));
//                for(int j=0; j<recDeteEB.rows(); j++) {
//                    if(recDeteEB.getShort("pindex",j)==loop && recDeteEB.getByte("detector",j)==17/*12*/) {
//                        int layer     = recDeteEB.getByte("layer",j);
//                        int paddle    = recDeteEB.getShort("component",j);
//                        double energy = recDeteEB.getFloat("energy",j);
////                System.out.println(ev + " " + pidCode + " " + recParticle.charge() + " " + recBankEB.getByte("charge", loop));
////                recBankEB.show();
//                        if(recParticle.charge()>0) {
//                            this.getDataGroup().getItem(1).getH1F("hi_pos_en_"+layer).fill(energy);
//                            this.getDataGroup().getItem(1).getH2F("hi_pos_en_p_"+layer).fill(recParticle.p(),energy);
//                            this.getDataGroup().getItem(1).getH2F("hi_pos_en_paddle_"+layer).fill(paddle*1.,energy);
//                        }
//                        else {
//                            this.getDataGroup().getItem(2).getH1F("hi_neg_en_"+layer).fill(energy);
//                            this.getDataGroup().getItem(2).getH2F("hi_neg_en_p_"+layer).fill(recParticle.p(),energy);
//                            this.getDataGroup().getItem(2).getH2F("hi_neg_en_paddle_"+layer).fill(paddle*1.,energy);
//                        }
//                    }
//                }
//            }
//        }
        if(recFtofHits!=null) {
            int nrows = recFtofHits.rows();
            for(int loop = 0; loop < nrows; loop++){
                int sector    = recFtofHits.getByte("sector", loop);
                int layer     = recFtofHits.getByte("layer", loop);
                int paddle    = recFtofHits.getShort("component", loop);
                int trk_id    = recFtofHits.getShort("trackid", loop);
                double energy = recFtofHits.getFloat("energy", loop);
                double time   = recFtofHits.getFloat("time",loop);
           	double tx     = recFtofHits.getFloat("tx", loop);
                double ty     = recFtofHits.getFloat("ty", loop);
                double tz     = recFtofHits.getFloat("tz", loop);
                double dx     = recFtofHits.getFloat("pathLengthThruBar", loop);
        	int adcId1    = recFtofHits.getShort("adc_idx1", loop);
                int adcId2    = recFtofHits.getShort("adc_idx2", loop);
                int tdcId1    = recFtofHits.getShort("tdc_idx1", loop);
                int tdcId2    = recFtofHits.getShort("tdc_idx2", loop);
                double adc1   = (double) ftofADC.getInt("ADC",adcId1);
                double adc2   = (double) ftofADC.getInt("ADC",adcId2);
                double adct1  = (double) ftofADC.getFloat("time",adcId1);
                double adct2  = (double) ftofADC.getFloat("time",adcId2);
                double tdc1   = (double) ftofTDC.getInt("TDC",tdcId1);
                double tdc2   = (double) ftofTDC.getInt("TDC",tdcId2);
                Vector3D trk = new Vector3D(tx,ty,tz);
                double angle = Math.toRadians(90-60*(sector-1));
        	trk.rotateZ(angle);
        	angle = Math.toRadians(25);
        	trk.rotateX(angle);
                this.getDataGroup().getItem(3).getH2F("hi_en_paddle_"+layer).fill(paddle*1.,energy);
                if(trk_id!=-1 && energy>1.5 && recHBTTrack!=null && recRunRF!=null) {
        	    int    q    = recHBTTrack.getByte("q",trk_id-1);
                    double c3x  = recHBTTrack.getFloat("c3_x",trk_id-1);
                    double c3y  = recHBTTrack.getFloat("c3_y",trk_id-1);
                    double c3z  = recHBTTrack.getFloat("c3_z",trk_id-1);
                    double p0x  = recHBTTrack.getFloat("p0_x",trk_id-1);
                    double p0y  = recHBTTrack.getFloat("p0_y",trk_id-1);
                    double p0z  = recHBTTrack.getFloat("p0_z",trk_id-1);
                    double vx   = recHBTTrack.getFloat("Vtx0_x",trk_id-1);
                    double vy   = recHBTTrack.getFloat("Vtx0_y",trk_id-1);
                    double vz   = recHBTTrack.getFloat("Vtx0_z",trk_id-1);
                    double chi2 = recHBTTrack.getFloat("chi2",trk_id-1)/recHBTTrack.getShort("ndf",trk_id-1);
//                    System.out.println(chi2);
                    Particle recParticle = new Particle(211,p0x,p0y,p0z,0,0,0);
                    double beta = recParticle.p()/Math.sqrt(recParticle.p()*recParticle.p()+recParticle.mass2());
                    double path = recHBTTrack.getFloat("pathlength",trk_id-1) + Math.sqrt((tx-c3x)*(tx-c3x)+(ty-c3y)*(ty-c3y)+(tz-c3z)*(tz-c3z));
                    double trf = 0;
                    for(int k = 0; k < recRunRF.rows(); k++){
                        if(recRunRF.getInt("id", k)==1) trf = recRunRF.getFloat("time",k);
                    }
                    double dt = (time - path/(beta*29.97) - trf + 120.5*2.004)%2.004-1.002;
                    if(/*Math.abs(trk.x())<10. && Math.abs(vz)<15.0 && adc1>1200 && adc2>1200 && */q==-1 && recParticle.p()>1) this.getDataGroup().getItem(4).getH2F("hi_rf_neg_paddle_"+layer).fill(paddle*1.,dt);
                    if(/*Math.abs(trk.x())<20. && Math.abs(vz)<15.0 && adc1>1200 && adc2>1200 && */q==1  && recParticle.p()>1) this.getDataGroup().getItem(5).getH2F("hi_rf_pos_paddle_"+layer).fill(paddle*1.,dt);
                    if(recEvenEB!=null) {
                        double startTime = recEvenEB.getFloat("STTime", 0);
                        if(startTime>-100) {
                            double betaTof = path/(time-startTime)/29.97;
                            double mass2   = Math.pow(recParticle.p()/betaTof, 2)-recParticle.p()*recParticle.p();
                            if(q==1)  {
                                this.getDataGroup().getItem(1).getH1F("hi_pos_en_"+layer).fill(energy/dx);
                                this.getDataGroup().getItem(1).getH2F("hi_pos_en_p_"+layer).fill(recParticle.p(),energy/dx);
                                this.getDataGroup().getItem(1).getH2F("hi_pos_en_paddle_"+layer).fill(paddle*1.,energy/dx);
                                this.getDataGroup().getItem(6).getH2F("hi_beta_pos_"+layer).fill(recParticle.p(),betaTof);
                                this.getDataGroup().getItem(7).getH1F("hi_mass_pos_"+layer).fill(mass2);
                            }
                            if(q==-1) {
                                this.getDataGroup().getItem(2).getH1F("hi_neg_en_"+layer).fill(energy/dx);
                                this.getDataGroup().getItem(2).getH2F("hi_neg_en_p_"+layer).fill(recParticle.p(),energy/dx);
                                this.getDataGroup().getItem(2).getH2F("hi_neg_en_paddle_"+layer).fill(paddle*1.,energy/dx);
                                this.getDataGroup().getItem(6).getH2F("hi_beta_neg_"+layer).fill(recParticle.p(),betaTof);
                                this.getDataGroup().getItem(7).getH1F("hi_mass_neg_"+layer).fill(mass2);
                            }
                        }
                    }
                }
            }
        }
        if(recFtofRaws!=null) {
            int nrows = recFtofRaws.rows();
            for(int loop = 0; loop < nrows; loop++){
                int layer    = recFtofRaws.getByte("layer", loop);
                int paddle   = recFtofRaws.getShort("component", loop);
                float tleft  = recFtofRaws.getFloat("time_left", loop);
                float tright = recFtofRaws.getFloat("time_right", loop);
                this.getDataGroup().getItem(3).getH2F("hi_time_paddle_"+layer).fill(tleft-tright,paddle*1.);
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
            H2F         h2    = this.getDataGroup().getItem(4).getH2F("hi_rf_neg_paddle_" + layer);
            GraphErrors meanX = this.getDataGroup().getItem(4).getGraph("g_rf_neg_paddle_" + layer);
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
            H2F         h2    = this.getDataGroup().getItem(5).getH2F("hi_rf_pos_paddle_" + layer);
            GraphErrors meanX = this.getDataGroup().getItem(5).getGraph("g_rf_pos_paddle_" + layer);
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
