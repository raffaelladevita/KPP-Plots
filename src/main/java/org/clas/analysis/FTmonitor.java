/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.analysis;

import java.util.ArrayList;
import org.clas.viewer.AnalysisMonitor;
import org.jlab.clas.physics.Particle;
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
public class FTmonitor extends AnalysisMonitor {
    
    double rfPeriod = 4.008;
    
    public FTmonitor(String name, ConstantsManager ccdb) {
        super(name,ccdb);
        this.setAnalysisTabNames("Hodoscope energy", "Hodoscope time", "Calorimeter", "Calorimeter time", "pi0");
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
   
        // Hodoscope
        DataGroup dc_hodo = new DataGroup(4,2);
        for(int layer=1; layer <= 2; layer++) {
            H1F hi_hodo_eall = new H1F("hi_hodo_eall_l"+layer, "E (MeV)", "Counts", 200, 0, 10); 
            hi_hodo_eall.setFillColor(4);
//            hi_hodo_eall.setOptStat("1111111");
            H1F hi_hodo_ematch = new H1F("hi_hodo_ematch_l"+layer,  "E (MeV)", "Counts", 200, 0, 10); 
            hi_hodo_ematch.setFillColor(3);
            F1D f_charge_landau = new F1D("Landau_" + layer,"[amp]*landau(x,[mean],[sigma])+[p0]+[p1]*x", 0.5*layer, 10.0);
            f_charge_landau.setParameter(0,0.0);
            f_charge_landau.setParameter(1,0.0);
            f_charge_landau.setParameter(2,1.0);
            f_charge_landau.setParameter(3,0.0);
            f_charge_landau.setParameter(4,0.0);  
            f_charge_landau.setOptStat(1111111); 
            f_charge_landau.setLineWidth(2);  
            H2F hi_hodo_ematch_2D = new H2F("hi_hodo_ematch_2D_l"+layer,  "hi_hodo_ematch_2D_l"+layer, 100, 0, 10 , 118, 0, 118); 
            hi_hodo_ematch_2D.setTitleX("E (MeV)"); 
            hi_hodo_ematch_2D.setTitleY("Tile");
            H1F hi_hodo_tmatch = new H1F("hi_hodo_tmatch_l"+layer,  "T-T_start (ns)", "Counts", 100, -50, 50); 
            hi_hodo_tmatch.setFillColor(3);
            H2F hi_hodo_tmatch_2D = new H2F("hi_hodo_tmatch_2D_l"+layer,  "hi_hodo_tmatch_2DD_l"+layer, 100, -50, 50, 118, 0, 118); 
            hi_hodo_tmatch_2D.setTitleX("E (MeV)"); 
            hi_hodo_tmatch_2D.setTitleY("Tile");
            dc_hodo.addDataSet(hi_hodo_eall,      0+4*(layer-1));
            dc_hodo.addDataSet(hi_hodo_ematch,    0+4*(layer-1));
            dc_hodo.addDataSet(f_charge_landau,   0+4*(layer-1));
            dc_hodo.addDataSet(hi_hodo_ematch_2D, 1+4*(layer-1));
            dc_hodo.addDataSet(hi_hodo_tmatch,    2+4*(layer-1));
            dc_hodo.addDataSet(hi_hodo_tmatch_2D, 3+4*(layer-1));
        }
        this.getDataGroup().add(dc_hodo, 0);
        // Calorimeter
        DataGroup dc_calo = new DataGroup(6,2);
        H1F hi_cal_nclusters = new H1F("hi_cal_nclusters", "N. Clusters", "Counts", 5, 0, 5);    
        hi_cal_nclusters.setFillColor(44);
        H1F hi_cal_clsize = new H1F("hi_cal_clsize", "Cluster Size", "Counts", 25, 0, 25); 
        hi_cal_clsize.setFillColor(44);
        H1F hi_cal_clsize_ch = new H1F("hi_cal_clsize_ch", "Cluster Size", "Counts", 25, 0, 25);  
        hi_cal_clsize_ch.setFillColor(44);
        H2F hi_cal_clsize_en = new H2F("hi_cal_clsize_en", " ", 25, 0, 25, 100, 0, 12); 
        hi_cal_clsize_en.setTitleX("Cluster size"); 
        hi_cal_clsize_en.setTitleY("E (GeV)");
        H1F hi_cal_e_all = new H1F("hi_cal_e_all", "E (GeV)", "Counts", 200, 0, 12); 
        hi_cal_e_all.setFillColor(4);
        H1F hi_cal_e_ch = new H1F("hi_cal_e_ch", "E (GeV)", "Counts", 200, 0, 12); 
        hi_cal_e_ch.setFillColor(2);
        H1F hi_cal_e_neu = new H1F("hi_cal_e_neu", "E (GeV)", "Counts", 200, 0, 12); 
        hi_cal_e_neu.setFillColor(3);
        H1F hi_cal_theta_ch = new H1F("hi_cal_theta_ch","#theta (deg)", "Counts", 100, 2,  6); 
        hi_cal_theta_ch.setFillColor(2);
        H1F hi_cal_phi_ch = new H1F("hi_cal_phi_ch", "#phi (deg)", "Counts", 100, -180,180); 
        hi_cal_phi_ch.setFillColor(2);
        H2F hi_cal_phi_e_ch = new H2F("hi_cal_phi_e_ch", "hi_cal_phi_e_ch", 100, -180,180, 200, 0, 12); 
        hi_cal_phi_e_ch.setTitleX("#phi (deg)");
        hi_cal_phi_e_ch.setTitleY("E (GeV)");
        H1F hi_cal_time_ch = new H1F("hi_cal_time_ch", "T-T_RF(ns)", "Counts", 100, -rfPeriod/2,rfPeriod/2); 
        hi_cal_time_ch.setFillColor(33);
        H1F hi_cal_time_cut_ch = new H1F("hi_cal_time_cut_ch", "T-T_RF(ns)", "Counts", 100, -rfPeriod/2,rfPeriod/2); 
        hi_cal_time_cut_ch.setFillColor(3);
        F1D ftime_ch = new F1D("ftime_ch", "[amp]*gaus(x,[mean],[sigma])", -1., 1.);
        ftime_ch.setParameter(0, 0.0);
        ftime_ch.setParameter(1, 0.0);
        ftime_ch.setParameter(2, 2.0);
        ftime_ch.setLineWidth(2);
        ftime_ch.setOptStat("1111");
        H2F hi_cal_time_e_ch = new H2F("hi_cal_time_e_ch", "hi_cal_time_e_ch", 100, 0., 12., 100, -rfPeriod/2,rfPeriod/2); 
        hi_cal_time_e_ch.setTitleX("E (GeV)");
        hi_cal_time_e_ch.setTitleY("T-T_RF (ns)");
        H2F hi_cal_time_theta_ch = new H2F("hi_cal_time_theta_ch", "hi_cal_time_theta_ch", 100, 2., 6., 100, -rfPeriod/2,rfPeriod/2); 
        hi_cal_time_theta_ch.setTitleX("#theta (deg)");
        hi_cal_time_theta_ch.setTitleY("T-T_RF (ns)");
        H1F hi_cal_time_neu = new H1F("hi_cal_time_neu", "T-T_start(ns)", "Counts", 100, -2,2); 
        hi_cal_time_neu.setFillColor(44);
        H1F hi_cal_time_cut_neu = new H1F("hi_cal_time_cut_neu", "T-T_start(ns)", "Counts", 100, -2,2); 
        hi_cal_time_cut_neu.setFillColor(4);
        F1D ftime_neu = new F1D("ftime_neu", "[amp]*gaus(x,[mean],[sigma])", -1., 1.);
        ftime_neu.setParameter(0, 0.0);
        ftime_neu.setParameter(1, 0.0);
        ftime_neu.setParameter(2, 2.0);
        ftime_neu.setLineWidth(2);
        ftime_neu.setOptStat("1111");
        H2F hi_cal_time_e_neu = new H2F("hi_cal_time_e_neu", "hi_cal_time_e_neu", 100, 0., 12., 100, -2,2);
        hi_cal_time_e_neu.setTitleX("E (GeV)");
        hi_cal_time_e_neu.setTitleY("T-T_start (ns)");
        H2F hi_cal_time_theta_neu = new H2F("hi_cal_time_theta_neu", "hi_cal_time_theta_neu", 100, 2., 6., 100, -2,2);
        hi_cal_time_theta_neu.setTitleX("#theta (deg)");
        hi_cal_time_theta_neu.setTitleY("T-T_start (ns)");
        dc_calo.addDataSet(hi_cal_nclusters, 0);
        dc_calo.addDataSet(hi_cal_clsize,    1);
        dc_calo.addDataSet(hi_cal_clsize_ch, 1);
        dc_calo.addDataSet(hi_cal_clsize_en, 2);
        dc_calo.addDataSet(hi_cal_e_all,     3);
        dc_calo.addDataSet(hi_cal_e_ch,      3);
        dc_calo.addDataSet(hi_cal_e_neu,     3);
        dc_calo.addDataSet(hi_cal_theta_ch,  4);
        dc_calo.addDataSet(hi_cal_phi_ch,    5);
        dc_calo.addDataSet(hi_cal_phi_e_ch,  5);
        dc_calo.addDataSet(hi_cal_time_ch,   6);
        dc_calo.addDataSet(hi_cal_time_cut_ch,   6);
        dc_calo.addDataSet(ftime_ch,         6);
        dc_calo.addDataSet(hi_cal_time_neu,  7);
        dc_calo.addDataSet(hi_cal_time_cut_neu,  7);
        dc_calo.addDataSet(ftime_neu,        7);
        dc_calo.addDataSet(hi_cal_time_e_ch, 8);
        dc_calo.addDataSet(hi_cal_time_e_neu,9);
        dc_calo.addDataSet(hi_cal_time_theta_ch, 10);
        dc_calo.addDataSet(hi_cal_time_theta_neu,11);
        this.getDataGroup().add(dc_calo, 1);
        // pi0
        H1F hpi0sum = new H1F("hpi0sum", 200,50., 250.);
        hpi0sum.setTitleX("M (MeV)");
        hpi0sum.setTitleY("Counts");
        hpi0sum.setTitle("2#gamma invariant mass");
        hpi0sum.setFillColor(3);
        F1D fpi0 = new F1D("fpi0", "[amp]*gaus(x,[mean],[sigma])+[p0]+[p1]*x", 80.,200.);
        fpi0.setParameter(0, 0.0);
        fpi0.setParameter(1, 140.0);
        fpi0.setParameter(2, 2.0);
        fpi0.setParameter(3, 0.0);
        fpi0.setParameter(4, 0.0);
        fpi0.setLineWidth(2);
        fpi0.setOptStat("1111111");
        H2F hmassangle = new H2F("hmassangle", 100, 0., 300., 100, 0., 6.);
        hmassangle.setTitleX("M (MeV)");
        hmassangle.setTitleY("Angle (deg)");
        hmassangle.setTitle("Angle vs. Mass");
        DataGroup dc_pi0 = new DataGroup(2,1);
        dc_pi0.addDataSet(hpi0sum,    0);
        dc_pi0.addDataSet(fpi0,       0);
        dc_pi0.addDataSet(hmassangle, 1);
        this.getDataGroup().add(dc_pi0, 2);
    }
    
    @Override
    public void plotHistos() {

        this.getAnalysisCanvas().getCanvas("Hodoscope energy").divide(2, 2);
        this.getAnalysisCanvas().getCanvas("Hodoscope energy").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Hodoscope energy").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Hodoscope time").divide(2, 2);
        this.getAnalysisCanvas().getCanvas("Hodoscope time").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Hodoscope time").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Calorimeter").divide(3, 2);
        this.getAnalysisCanvas().getCanvas("Calorimeter").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Calorimeter").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Calorimeter time").divide(3, 2);
        this.getAnalysisCanvas().getCanvas("Calorimeter time").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Calorimeter time").setGridY(false);
        this.getAnalysisCanvas().getCanvas("pi0").divide(2,1);
        this.getAnalysisCanvas().getCanvas("pi0").setGridX(false);
        this.getAnalysisCanvas().getCanvas("pi0").setGridY(false);
        
        for(int layer=1; layer <= 2; layer++) {
            this.getAnalysisCanvas().getCanvas("Hodoscope energy").cd(0+2*(layer-1));
//        this.getAnalysisCanvas().getCanvas("Beta").getPad(0).getAxisZ().setLog(true);
            this.getAnalysisCanvas().getCanvas("Hodoscope energy").draw(this.getDataGroup().getItem(0).getH1F("hi_hodo_eall_l"+layer));
            this.getAnalysisCanvas().getCanvas("Hodoscope energy").draw(this.getDataGroup().getItem(0).getH1F("hi_hodo_ematch_l"+layer),"same");
            this.getAnalysisCanvas().getCanvas("Hodoscope energy").draw(this.getDataGroup().getItem(0).getF1D("Landau_"+layer),"same");
            this.getAnalysisCanvas().getCanvas("Hodoscope energy").cd(1+2*(layer-1));
            this.getAnalysisCanvas().getCanvas("Hodoscope energy").getPad(1+2*(layer-1)).getAxisZ().setLog(true);
            this.getAnalysisCanvas().getCanvas("Hodoscope energy").draw(this.getDataGroup().getItem(0).getH2F("hi_hodo_ematch_2D_l"+layer));
            this.getAnalysisCanvas().getCanvas("Hodoscope time").cd(0+2*(layer-1));
//        this.getAnalysisCanvas().getCanvas("Beta").getPad(0).getAxisZ().setLog(true);
            this.getAnalysisCanvas().getCanvas("Hodoscope time").draw(this.getDataGroup().getItem(0).getH1F("hi_hodo_tmatch_l"+layer));
            this.getAnalysisCanvas().getCanvas("Hodoscope time").cd(1+2*(layer-1));
            this.getAnalysisCanvas().getCanvas("Hodoscope time").getPad(1+2*(layer-1)).getAxisZ().setLog(true);
            this.getAnalysisCanvas().getCanvas("Hodoscope time").draw(this.getDataGroup().getItem(0).getH2F("hi_hodo_tmatch_2D_l"+layer));
        }
        this.getAnalysisCanvas().getCanvas("Calorimeter").cd(0);
        this.getAnalysisCanvas().getCanvas("Calorimeter").getPad(0).getAxisY().setLog(true);
        this.getAnalysisCanvas().getCanvas("Calorimeter").draw(this.getDataGroup().getItem(1).getH1F("hi_cal_nclusters"));
        this.getAnalysisCanvas().getCanvas("Calorimeter").cd(1);
        this.getAnalysisCanvas().getCanvas("Calorimeter").getPad(1).getAxisY().setLog(true);
        this.getAnalysisCanvas().getCanvas("Calorimeter").draw(this.getDataGroup().getItem(1).getH1F("hi_cal_clsize"));
        this.getAnalysisCanvas().getCanvas("Calorimeter").draw(this.getDataGroup().getItem(1).getH1F("hi_cal_clsize_ch"),"same");
        this.getAnalysisCanvas().getCanvas("Calorimeter").cd(2);
        this.getAnalysisCanvas().getCanvas("Calorimeter").getPad(2).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("Calorimeter").draw(this.getDataGroup().getItem(1).getH2F("hi_cal_clsize_en"));
        this.getAnalysisCanvas().getCanvas("Calorimeter").cd(3);
        this.getAnalysisCanvas().getCanvas("Calorimeter").draw(this.getDataGroup().getItem(1).getH1F("hi_cal_e_all"));
        this.getAnalysisCanvas().getCanvas("Calorimeter").draw(this.getDataGroup().getItem(1).getH1F("hi_cal_e_ch"),"same");
        this.getAnalysisCanvas().getCanvas("Calorimeter").draw(this.getDataGroup().getItem(1).getH1F("hi_cal_e_neu"),"same");
        this.getAnalysisCanvas().getCanvas("Calorimeter").cd(4);
        this.getAnalysisCanvas().getCanvas("Calorimeter").draw(this.getDataGroup().getItem(1).getH1F("hi_cal_theta_ch"));
        this.getAnalysisCanvas().getCanvas("Calorimeter").cd(5);
        this.getAnalysisCanvas().getCanvas("Calorimeter").draw(this.getDataGroup().getItem(1).getH1F("hi_cal_phi_ch"));
//        this.getAnalysisCanvas().getCanvas("Calorimeter").draw(this.getDataGroup().getItem(1).getH2F("hi_cal_phi_e_ch"));
        this.getAnalysisCanvas().getCanvas("Calorimeter time").cd(0);
        this.getAnalysisCanvas().getCanvas("Calorimeter time").draw(this.getDataGroup().getItem(1).getH1F("hi_cal_time_ch"));
        this.getAnalysisCanvas().getCanvas("Calorimeter time").draw(this.getDataGroup().getItem(1).getH1F("hi_cal_time_cut_ch"),"same");
        this.getAnalysisCanvas().getCanvas("Calorimeter time").draw(this.getDataGroup().getItem(1).getF1D("ftime_ch"),"same");
        this.getAnalysisCanvas().getCanvas("Calorimeter time").cd(1);
        this.getAnalysisCanvas().getCanvas("Calorimeter time").draw(this.getDataGroup().getItem(1).getH2F("hi_cal_time_e_ch"));
        this.getAnalysisCanvas().getCanvas("Calorimeter time").cd(2);
        this.getAnalysisCanvas().getCanvas("Calorimeter time").draw(this.getDataGroup().getItem(1).getH2F("hi_cal_time_theta_ch"));
        this.getAnalysisCanvas().getCanvas("Calorimeter time").cd(3);
        this.getAnalysisCanvas().getCanvas("Calorimeter time").draw(this.getDataGroup().getItem(1).getH1F("hi_cal_time_neu"));
        this.getAnalysisCanvas().getCanvas("Calorimeter time").draw(this.getDataGroup().getItem(1).getH1F("hi_cal_time_cut_neu"),"same");
        this.getAnalysisCanvas().getCanvas("Calorimeter time").draw(this.getDataGroup().getItem(1).getF1D("ftime_neu"),"same");
        this.getAnalysisCanvas().getCanvas("Calorimeter time").cd(4);
        this.getAnalysisCanvas().getCanvas("Calorimeter time").draw(this.getDataGroup().getItem(1).getH2F("hi_cal_time_e_neu"));
        this.getAnalysisCanvas().getCanvas("Calorimeter time").cd(5);
        this.getAnalysisCanvas().getCanvas("Calorimeter time").draw(this.getDataGroup().getItem(1).getH2F("hi_cal_time_theta_neu"));
        
        this.getAnalysisCanvas().getCanvas("pi0").cd(0);
        this.getAnalysisCanvas().getCanvas("pi0").draw(this.getDataGroup().getItem(2).getH1F("hpi0sum"));
        this.getAnalysisCanvas().getCanvas("pi0").draw(this.getDataGroup().getItem(2).getF1D("fpi0"),"same");
        this.getAnalysisCanvas().getCanvas("pi0").cd(1);
        this.getAnalysisCanvas().getCanvas("pi0").draw(this.getDataGroup().getItem(2).getH2F("hmassangle"));
    }
        
    @Override
    public void processEvent(DataEvent event) {
        // process event info and save into data group        
        // event builder
        DataBank recRun    = null;
        DataBank recBankEB = null;
        DataBank recEvenEB = null;
        DataBank ftParticles = null;
        DataBank ftCalClusters = null;
        DataBank ftCalHits = null;
        DataBank ftHodoClusters = null;
        DataBank ftHodoHits = null;
        if(event.hasBank("RUN::config"))            recRun      = event.getBank("RUN::config");
        if(event.hasBank("REC::Particle"))          recBankEB   = event.getBank("REC::Particle");
        if(event.hasBank("REC::Event"))             recEvenEB   = event.getBank("REC::Event");
        if(event.hasBank("FT::particles"))          ftParticles = event.getBank("FT::particles");
        if(event.hasBank("FTCAL::clusters"))      ftCalClusters = event.getBank("FTCAL::clusters");
        if(event.hasBank("FTCAL::hits"))              ftCalHits = event.getBank("FTCAL::hits");
        if(event.hasBank("FTHODO::clusters"))    ftHodoClusters = event.getBank("FTHODO::clusters");
        if(event.hasBank("FTHODO::hits"))            ftHodoHits = event.getBank("FTHODO::hits");
        if(recRun == null) return;
        int ev  = recRun.getInt("event",0);
        int run = recRun.getInt("run",0);
        IndexedTable rfConfig = this.getCcdb().getConstants(run, "/calibration/eb/rf/config");
        if(this.rfPeriod!=rfConfig.getDoubleValue("clock", 1,1,1)) {
            this.rfPeriod = rfConfig.getDoubleValue("clock", 1,1,1);
            this.resetEventListener();
        }
        // Decoding Trigger Bits
        boolean[] trigger_bits = new boolean[32];
        if (event.hasBank("RUN::config")) {
            DataBank bank = event.getBank("RUN::config");
            long TriggerWord = bank.getLong("trigger", 0) & 0xFFFFFFFF;
            for (int i = 31; i >= 0; i--) {
                trigger_bits[i] = (TriggerWord & (1 << i)) != 0;
            }
        }
        // get event start time
        double startTime=-1000;
        double rfTime=-1000;
        if(recEvenEB!=null) {
            startTime = recEvenEB.getFloat("startTime", 0);
            rfTime    = recEvenEB.getFloat("RFTime", 0);
        }
        // get trigger particle
        int trigger=0;
        if(recBankEB!=null) trigger = recBankEB.getInt("pid", 0);
     
        if (ftParticles != null /*&& trigger_bits[25]*/) {
            Particle p1 = null;
            Particle p2 = null;
            ArrayList<Particle> gammas = new ArrayList();
            this.getDataGroup().getItem(1).getH1F("hi_cal_nclusters").fill(ftParticles.rows());
            for (int loop = 0; loop < ftParticles.rows(); loop++) {
                int    charge = ftParticles.getByte("charge", loop);
                double energy = ftParticles.getFloat("energy", loop);
                double   time = ftParticles.getFloat("time", loop);
                double     cx = ftParticles.getFloat("cx", loop);
                double     cy = ftParticles.getFloat("cy", loop);
                double     cz = ftParticles.getFloat("cz", loop);
                int     calID = ftParticles.getShort("calID", loop);
                int    hodoID = ftParticles.getShort("hodoID", loop);
                
	        double energyR  = 0; 
	        int    size     = 0;
                double path     = 0;
	        for(int i=0; i<ftCalClusters.rows(); i++) {
	            if(calID == ftCalClusters.getShort("id", i)) {
                        energyR  = ftCalClusters.getFloat("recEnergy", i);
                        size     = ftCalClusters.getInt("size", i);
                        double x = ftCalClusters.getFloat("x", i);
                        double y = ftCalClusters.getFloat("y", i);
                        double z = ftCalClusters.getFloat("z", i);
                        path     = Math.sqrt(x*x+y*y+z*z);
                        time     = ftCalClusters.getFloat("time", i)-path/29.97;
		    }
	        }
//	        for(int i=0; i<ftCalHits.rows(); i++) {
//                        energyR  = ftCalHits.getFloat("energy", i);
//                        double x = ftCalHits.getFloat("x", i);
//                        double y = ftCalHits.getFloat("y", i);
//                        double z = ftCalHits.getFloat("z", i);
//                        path     = Math.sqrt(x*x+y*y+z*z);
//                        time     = ftCalHits.getFloat("time", i)-path/29.97;
//                        this.getDataGroup().getItem(1).getH1F("hi_cal_time_ch").fill(time-startTime);
//                        if(energy>2) this.getDataGroup().getItem(1).getH1F("hi_cal_time_cut_ch").fill(time-startTime);
//                        this.getDataGroup().getItem(1).getH2F("hi_cal_time_e_ch").fill(energy,(time-startTime));
//	        }
                
                this.getDataGroup().getItem(1).getH1F("hi_cal_clsize").fill(size);
                this.getDataGroup().getItem(1).getH1F("hi_cal_e_all").fill(energy);
                this.getDataGroup().getItem(1).getH2F("hi_cal_clsize_en").fill(size, energy);

                if (charge != 0) {
                    this.getDataGroup().getItem(1).getH1F("hi_cal_clsize_ch").fill(size);
                    this.getDataGroup().getItem(1).getH1F("hi_cal_e_ch").fill(energy);
                    this.getDataGroup().getItem(1).getH1F("hi_cal_theta_ch").fill(Math.toDegrees(Math.acos(cz)));
                    if(energy>0.5)this.getDataGroup().getItem(1).getH1F("hi_cal_phi_ch").fill(Math.toDegrees(Math.atan2(cy,cx)));
                    this.getDataGroup().getItem(1).getH2F("hi_cal_phi_e_ch").fill(Math.toDegrees(Math.atan2(cy,cx)),energy);
                    if(rfTime!=-1000) {
                        if(energy>0.5) this.getDataGroup().getItem(1).getH1F("hi_cal_time_ch").fill((time-rfTime+1000.5*rfPeriod)%rfPeriod-0.5*rfPeriod);
                        if(energy>2.0) this.getDataGroup().getItem(1).getH1F("hi_cal_time_cut_ch").fill((time-rfTime+1000.5*rfPeriod)%rfPeriod-0.5*rfPeriod);
                        this.getDataGroup().getItem(1).getH2F("hi_cal_time_e_ch").fill(energy,(time-rfTime+1000.5*rfPeriod)%rfPeriod-0.5*rfPeriod);
                        this.getDataGroup().getItem(1).getH2F("hi_cal_time_theta_ch").fill(Math.toDegrees(Math.acos(cz)),(time-rfTime+1000.5*rfPeriod)%rfPeriod-0.5*rfPeriod);
                    }
                }
                else {
                    Particle recParticle = new Particle(22, energy*cx, energy*cy, energy*cz, 0,0,0);
                    gammas.add(recParticle);
                    this.getDataGroup().getItem(1).getH1F("hi_cal_e_neu").fill(energy); 
                    if(startTime!=-1000 && trigger==11) {
                        this.getDataGroup().getItem(1).getH1F("hi_cal_time_neu").fill(time-startTime);
                        if(energy>2) this.getDataGroup().getItem(1).getH1F("hi_cal_time_cut_neu").fill(time-startTime);
                        this.getDataGroup().getItem(1).getH2F("hi_cal_time_e_neu").fill(energy,time-startTime);
                        this.getDataGroup().getItem(1).getH2F("hi_cal_time_theta_neu").fill(Math.toDegrees(Math.acos(cz)),time-startTime);
                    }
                }
            }
            if(gammas.size()>=2) {
                for (int i1 = 0; i1 < gammas.size(); i1++) {
                    for (int i2 = i1 + 1; i2 < gammas.size(); i2++) {
                        Particle partGamma1 = gammas.get(i1);
                        Particle partGamma2 = gammas.get(i2);
                        Particle partPi0 = new Particle();
                        partPi0.copy(partGamma1);
                        partPi0.combine(partGamma2, +1);
                        double invmass = Math.sqrt(partPi0.mass2());
                        double x = (partGamma1.p() - partGamma2.p()) / (partGamma1.p() + partGamma2.p());
                        double angle = Math.toDegrees(Math.acos(partGamma1.cosTheta(partGamma2)));
                        if(angle>1.5) this.getDataGroup().getItem(2).getH1F("hpi0sum").fill(invmass*1000);
                        this.getDataGroup().getItem(2).getH2F("hmassangle").fill(invmass*1000, angle);
                    }
                }
            }       
            if(ftHodoHits!= null) {
               for(int i=0; i<ftHodoHits.rows(); i++) {
                   int hodoC = ftHodoHits.getShort("clusterID",i);
                   int hodoS = ftHodoHits.getByte("sector",i);
                   int hodoL = ftHodoHits.getByte("layer",i);
                   int component = ftHodoHits.getShort("component",i);
                   int tile = -1;
                   switch (hodoS) {
                       case 1:
                          tile = component + 0;
                          break;
                       case 2:
                          tile = component + 9;
                          break;
                       case 3:
                          tile = component + 29;
                          break;
                       case 4:
                          tile = component + 38;
                          break;
                       case 5:
                          tile = component + 58;
                          break;
                       case 6:
                          tile = component + 67;
                          break;
                       case 7:
                          tile = component + 87;
                          break;
                       case 8:
                          tile = component + 96;
                          break;
                       default:
                          tile = -1;
                          break;
                    }
                    double hodoHitE = ftHodoHits.getFloat("energy",i);
                    double hodoHitT = ftHodoHits.getFloat("time",i);
                    double hodoHitX = ftHodoHits.getFloat("x",i);
                    double hodoHitY = ftHodoHits.getFloat("y",i);
                    double hodoHitZ = ftHodoHits.getFloat("z",i);
                    double path = Math.sqrt(hodoHitX*hodoHitX+hodoHitY*hodoHitY+hodoHitZ*hodoHitZ);
                    int   clusterId = ftHodoHits.getShort("clusterID",i);
                    this.getDataGroup().getItem(0).getH1F("hi_hodo_eall_l"+hodoL).fill(hodoHitE);
//                    this.getDataGroup().getItem(0).getH2F("hi_hodo_eall_2D_l"+hodoL).fill(hodoHitE,tile);
                    for(int j=0; j<ftHodoClusters.rows(); j++) {
                        if(clusterId==ftHodoClusters.getShort("id", j) && ftHodoClusters.getShort("size", j)>1) {
                            this.getDataGroup().getItem(0).getH1F("hi_hodo_ematch_l"+hodoL).fill(hodoHitE);
                            this.getDataGroup().getItem(0).getH2F("hi_hodo_ematch_2D_l"+hodoL).fill(hodoHitE,tile);
                            if(startTime > -100) {
                                this.getDataGroup().getItem(0).getH1F("hi_hodo_tmatch_l"+hodoL).fill(hodoHitT-path/29.97-startTime);
                                this.getDataGroup().getItem(0).getH2F("hi_hodo_tmatch_2D_l"+hodoL).fill(hodoHitT-path/29.97-startTime,tile); 
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
        // fit hodoscope charge
        for(int layer=1; layer <= 2; layer++) {
            H1F h1 = this.getDataGroup().getItem(0).getH1F("hi_hodo_ematch_l"+layer);
            F1D f1 = this.getDataGroup().getItem(0).getF1D("Landau_"+layer);
            this.initLandauFitPar(h1, f1);
            DataFitter.fit(f1,h1,"LRQ");
            h1.setFunction(null); 
        }
        // fit calorimeter time
        H1F htime = this.getDataGroup().getItem(1).getH1F("hi_cal_time_cut_ch");
        F1D ftime = this.getDataGroup().getItem(1).getF1D("ftime_ch");
        this.initTimeGaussFitPar(ftime,htime);
        DataFitter.fit(ftime,htime,"LQ");
        htime.setFunction(null);
        htime = this.getDataGroup().getItem(1).getH1F("hi_cal_time_cut_neu");
        ftime = this.getDataGroup().getItem(1).getF1D("ftime_neu");
        this.initTimeGaussFitPar(ftime,htime);
        DataFitter.fit(ftime,htime,"LQ");
        htime.setFunction(null);
        // fit pi0 mass
        H1F h1p = this.getDataGroup().getItem(2).getH1F("hpi0sum");
        F1D f1p = this.getDataGroup().getItem(2).getF1D("fpi0");
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
    
    private void initLandauFitPar(H1F hcharge, F1D fcharge) {
        double hAmp  = hcharge.getBinContent(hcharge.getMaximumBin());
        double hMean = hcharge.getAxis().getBinCenter(hcharge.getMaximumBin());
        double hRMS  = hcharge.getRMS(); //ns
        fcharge.setRange(fcharge.getRange().getMin(), hMean*2.0);
        fcharge.setParameter(0, hAmp);
        fcharge.setParLimits(0, 0.5*hAmp, 1.5*hAmp); 
        fcharge.setParameter(1, hMean);
        fcharge.setParLimits(1, 0.8*hMean, 1.2*hMean);//Changed from 5-30        
        fcharge.setParameter(2, 0.3);//Changed from 2
        fcharge.setParLimits(2, 0.1, 1);//Changed from 0.5-10
        fcharge.setParameter(3, 0.2*hAmp);
//        fcharge.setParLimits(3,  0.0, 10000000.); 
        fcharge.setParameter(4, -0.3);//Changed from -0.2
//        fcharge.setParLimits(4,  0.0, 3.0); //Changed from -10-0
    }

    private void initTimeGaussFitPar(F1D ftime, H1F htime) {
        double hAmp  = htime.getBinContent(htime.getMaximumBin());
        double hMean = htime.getAxis().getBinCenter(htime.getMaximumBin());
        double hRMS  = htime.getRMS(); //ns
        double rangeMin = (hMean - (3*hRMS)); 
        double rangeMax = (hMean + (3*hRMS));  
        double pm = hRMS*3;
        ftime.setRange(rangeMin, rangeMax);
        ftime.setParameter(0, hAmp);
        ftime.setParLimits(0, hAmp*0.8, hAmp*1.2);
        ftime.setParameter(1, hMean);
        ftime.setParLimits(1, hMean-pm, hMean+(pm));
        ftime.setParameter(2, 0.2);
        ftime.setParLimits(2, 0.1*hRMS, 0.8*hRMS);
    }    

}
