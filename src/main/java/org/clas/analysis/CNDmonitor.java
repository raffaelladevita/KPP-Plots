/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.analysis;

import java.util.ArrayList;
import org.clas.viewer.AnalysisMonitor;
import org.jlab.clas.physics.Particle;
import org.jlab.detector.base.DetectorType;
import org.jlab.detector.calib.utils.ConstantsManager;
import org.jlab.groot.data.GraphErrors;
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
public class CNDmonitor extends AnalysisMonitor {
    
    private double rfPeriod = 4.008;
    private int    nPaddle = 144;
    private double tdcConv = 0.023456;
        
    public CNDmonitor(String name, ConstantsManager ccdb) {
        super(name,ccdb);
        this.setAnalysisTabNames("MIPs", "ADC-TDC", "dE/dx", "Residuals", "RF","Beta and Mass");
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
        H2F hi_time_paddle = new H2F("hi_time_paddle", "hi_time_paddle", nPaddle, 1, nPaddle+1., 200, -20., 20.);  
        hi_time_paddle.setTitleX("Counter"); 
        hi_time_paddle.setTitleY("Up-Down Time (ns)");
        dc_mips.addDataSet(hi_en_paddle,  0);
        dc_mips.addDataSet(hi_time_paddle,1);        
        this.getDataGroup().add(dc_mips, 2);  
        // Space residuals
        DataGroup dc_residuals = new DataGroup(3,2);
        H1F hi_z_hit = new H1F("hi_z_hit", "hi_z_hit", 100, -60., 40.);  
        hi_z_hit.setTitleX("CND hit z (cm)");
        hi_z_hit.setTitleY("Counts"); 
        H1F hi_z_track = new H1F("hi_z_track", "hi_z_track", 100, -60., 40.);  
        hi_z_track.setTitleX("CND track z (cm)");
        hi_z_track.setTitleY("Counts"); 
        H2F hi_z_hit_track = new H2F("hi_z_hit_track", "hi_z_hit_track", 100, -60., 40., 100, -60., 40.);  
        hi_z_hit_track.setTitleX("CND hit z (cm)");
        hi_z_hit_track.setTitleY("Track intersect z (cm)"); 
        H2F hi_z_residual = new H2F("hi_z_residual", "hi_z_residual", 200, -30., 30., nPaddle, 1, nPaddle+1.);  
        hi_z_residual.setTitleX("CND -Track z (cm)");
        hi_z_residual.setTitleY("Paddle"); 
        H2F hi_rz_hit = new H2F("hi_rz_hit", "hi_rz_hit", 100, -60., 40., 100, 20., 40.);  
        hi_rz_hit.setTitleX("CND hit z (cm)");
        hi_rz_hit.setTitleY("CND hit r (cm)"); 
        H2F hi_rz_track = new H2F("hi_rz_track", "hi_rz_track", 100, -60., 40., 100, 20., 40.);  
        hi_rz_track.setTitleX("Track z (cm)");
        hi_rz_track.setTitleY("Tracl r (cm)"); 
        dc_residuals.addDataSet(hi_z_hit,       0);
        dc_residuals.addDataSet(hi_rz_hit,      1);
        dc_residuals.addDataSet(hi_z_hit_track, 2);
        dc_residuals.addDataSet(hi_z_track,     3);
        dc_residuals.addDataSet(hi_rz_track,    4);
        dc_residuals.addDataSet(hi_z_residual,  5);        
        this.getDataGroup().add(dc_residuals, 3);   
        // RF offsets
        DataGroup dc_rf = new DataGroup(3,2);
        H2F hi_rf_neg_paddle = new H2F("hi_rf_neg_paddle", "hi_rf_neg_paddle", nPaddle, 1, nPaddle+1., 100, -this.rfPeriod/2, this.rfPeriod/2);  
        hi_rf_neg_paddle.setTitleX("Counter");
        hi_rf_neg_paddle.setTitleY("Vertex Time (ns)"); 
        H2F hi_rf_pos_paddle = new H2F("hi_rf_pos_paddle", "hi_rf_pos_paddle", nPaddle, 1, nPaddle+1., 100, -this.rfPeriod/2, this.rfPeriod/2);  
        hi_rf_pos_paddle.setTitleX("Counter");
        hi_rf_pos_paddle.setTitleY("Vertex Time (ns)"); 
        GraphErrors g_rf_neg_paddle = new GraphErrors("g_rf_neg_paddle");
        g_rf_neg_paddle.setTitleX("RF offset (ns)"); 
        g_rf_neg_paddle.setTitleY("Counter");
        GraphErrors g_rf_pos_paddle = new GraphErrors("g_rf_pos_paddle");
        g_rf_pos_paddle.setTitleX("RF offset (ns)"); 
        g_rf_pos_paddle.setTitleY("Counter");
        for(int i=1; i<=nPaddle; i++) g_rf_pos_paddle.addPoint((double) i, 0, 0, 0);
        H1F hi_rf_neg = new H1F("hi_rf_neg", "hi_rf_neg", 200, -this.rfPeriod/2, this.rfPeriod/2);  
        hi_rf_neg.setTitleX("Vertex Time (ns)");
        hi_rf_neg.setTitleY("Counts"); 
        F1D f1_rf_neg = new F1D("f1_rf_neg","[amp]*gaus(x,[mean],[sigma])+[p0]", -1.0, 1.0);
        f1_rf_neg.setParameter(0, 0);
        f1_rf_neg.setParameter(1, 0);
        f1_rf_neg.setParameter(2, 0.08);
        f1_rf_neg.setParameter(3, 0);
        f1_rf_neg.setLineWidth(2);
        f1_rf_neg.setLineColor(2);
        f1_rf_neg.setOptStat("1111");
        H1F hi_rf_pos = new H1F("hi_rf_pos", "hi_rf_pos", 200, -this.rfPeriod/2, this.rfPeriod/2);  
        hi_rf_pos.setTitleX("Vertex Time (ns)");
        hi_rf_pos.setTitleY("Counts"); 
        F1D f1_rf_pos = new F1D("f1_rf_pos","[amp]*gaus(x,[mean],[sigma])+[p0]", -1.0, 1.0);
        f1_rf_pos.setParameter(0, 0);
        f1_rf_pos.setParameter(1, 0);
        f1_rf_pos.setParameter(2, 0.08);
        f1_rf_pos.setParameter(3, 0);
        f1_rf_pos.setLineWidth(2);
        f1_rf_pos.setLineColor(2);
        f1_rf_pos.setOptStat("1111");
        dc_rf.addDataSet(hi_rf_pos_paddle, 0);
        dc_rf.addDataSet(hi_rf_neg_paddle, 1);
        dc_rf.addDataSet(hi_rf_pos,        2);
        dc_rf.addDataSet(f1_rf_pos,        2);
        dc_rf.addDataSet(g_rf_neg_paddle,  3);
        dc_rf.addDataSet(g_rf_pos_paddle,  4);        
        dc_rf.addDataSet(hi_rf_neg,        5);
        dc_rf.addDataSet(f1_rf_neg,        5);
        this.getDataGroup().add(dc_rf, 4);   
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
        F1D cpion = new F1D("cpion","x/sqrt(x*x+0.1396*0.1396)", 0.2, 2.0);
        F1D ckaon = new F1D("ckaon","x/sqrt(x*x+0.4935*0.4935)", 0.2, 2.0);
        F1D cprot = new F1D("cprot","x/sqrt(x*x+0.9383*0.9383)", 0.2, 2.0);
        F1D cdeut = new F1D("cdeut","x/sqrt(x*x+1.8756*1.8756)", 0.2, 2.0);
        dc_beta.addDataSet(hi_mass_pos,  0);
        dc_beta.addDataSet(hi_mass_neg,  1);
        dc_beta.addDataSet(hi_beta_pos,  2);
        dc_beta.addDataSet(cpion,        2);
        dc_beta.addDataSet(ckaon,        2);
        dc_beta.addDataSet(cprot,        2);
        dc_beta.addDataSet(cdeut,        2);  
        dc_beta.addDataSet(hi_beta_neg,  3);
        this.getDataGroup().add(dc_beta, 5);  
        // adc-tdc matching
        H2F hi_adctdc_left= new H2F("hi_adctdc_left", "hi_adctdc_left", 100, 0., 400, 100, 100., 500.);  
        hi_adctdc_left.setTitleX("ADC time (ns)"); 
        hi_adctdc_left.setTitleY("TDC time (ns)");
        H1F hi_adc_tdc_left= new H1F("hi_adc_tdc_left", "hi_adc_tdc_left", 100, 0., 200);  
        hi_adc_tdc_left.setTitleX("TDC - ADC time (ns)"); 
        hi_adc_tdc_left.setTitleY("Counts");
        H2F hi_adc_tdc_paddle_left= new H2F("hi_adc_tdc_paddle_left", "hi_adc_tdc_paddle_left", 100, 0., 200, nPaddle, 1, nPaddle+1.);  
        hi_adc_tdc_paddle_left.setTitleX("TDC - ADC time (ns)"); 
        hi_adc_tdc_paddle_left.setTitleY("Paddle");
        H2F hi_adctdc_right= new H2F("hi_adctdc_right", "hi_adctdc_right", 100, 0., 400, 100, 100., 500.);  
        hi_adctdc_right.setTitleX("ADC time (ns)"); 
        hi_adctdc_right.setTitleY("TDC time (ns)");
        H1F hi_adc_tdc_right= new H1F("hi_adc_tdc_right", "hi_adc_tdc_right", 100, 0., 200);  
        hi_adc_tdc_right.setTitleX("TDC - ADC time (ns)"); 
        hi_adc_tdc_right.setTitleY("Counts");
        H2F hi_adc_tdc_paddle_right= new H2F("hi_adc_tdc_paddle_right", "hi_adc_tdc_paddle_right", 100, 0., 200, nPaddle, 1, nPaddle+1.);  
        hi_adc_tdc_paddle_right.setTitleX("TDC - ADC time (ns)"); 
        hi_adc_tdc_paddle_right.setTitleY("Paddle");
        DataGroup dc_adctdc = new DataGroup(2,3);
        dc_adctdc.addDataSet(hi_adctdc_left,0);        
        dc_adctdc.addDataSet(hi_adc_tdc_left,1);        
        dc_adctdc.addDataSet(hi_adc_tdc_paddle_left,2);        
        dc_adctdc.addDataSet(hi_adctdc_right,3);        
        dc_adctdc.addDataSet(hi_adc_tdc_right,4);        
        dc_adctdc.addDataSet(hi_adc_tdc_paddle_right,5);        
        this.getDataGroup().add(dc_adctdc, 6);  
    }
    
    @Override
    public void plotHistos() {
        this.getAnalysisCanvas().getCanvas("MIPs").divide(2,1);
        this.getAnalysisCanvas().getCanvas("MIPs").setGridX(false);
        this.getAnalysisCanvas().getCanvas("MIPs").setGridY(false);
        this.getAnalysisCanvas().getCanvas("ADC-TDC").divide(3,2);
        this.getAnalysisCanvas().getCanvas("ADC-TDC").setGridX(false);
        this.getAnalysisCanvas().getCanvas("ADC-TDC").setGridY(false);
        this.getAnalysisCanvas().getCanvas("RF").divide(3, 2);
        this.getAnalysisCanvas().getCanvas("RF").setGridX(false);
        this.getAnalysisCanvas().getCanvas("RF").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Residuals").divide(3, 2);
        this.getAnalysisCanvas().getCanvas("Residuals").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Residuals").setGridY(false);
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
        this.getAnalysisCanvas().getCanvas("ADC-TDC").cd(0);
        this.getAnalysisCanvas().getCanvas("ADC-TDC").getPad(0).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("ADC-TDC").draw(this.getDataGroup().getItem(6).getH2F("hi_adctdc_left"));
        this.getAnalysisCanvas().getCanvas("ADC-TDC").cd(1);
        this.getAnalysisCanvas().getCanvas("ADC-TDC").getPad(1).getAxisY().setLog(true);
        this.getAnalysisCanvas().getCanvas("ADC-TDC").draw(this.getDataGroup().getItem(6).getH1F("hi_adc_tdc_left"));
        this.getAnalysisCanvas().getCanvas("ADC-TDC").cd(2);
        this.getAnalysisCanvas().getCanvas("ADC-TDC").getPad(2).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("ADC-TDC").draw(this.getDataGroup().getItem(6).getH2F("hi_adc_tdc_paddle_left"));
        this.getAnalysisCanvas().getCanvas("ADC-TDC").cd(3);
        this.getAnalysisCanvas().getCanvas("ADC-TDC").getPad(3).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("ADC-TDC").draw(this.getDataGroup().getItem(6).getH2F("hi_adctdc_right"));
        this.getAnalysisCanvas().getCanvas("ADC-TDC").cd(4);
        this.getAnalysisCanvas().getCanvas("ADC-TDC").getPad(4).getAxisY().setLog(true);
        this.getAnalysisCanvas().getCanvas("ADC-TDC").draw(this.getDataGroup().getItem(6).getH1F("hi_adc_tdc_right"));
        this.getAnalysisCanvas().getCanvas("ADC-TDC").cd(5);
        this.getAnalysisCanvas().getCanvas("ADC-TDC").getPad(5).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("ADC-TDC").draw(this.getDataGroup().getItem(6).getH2F("hi_adc_tdc_paddle_right"));
        this.getAnalysisCanvas().getCanvas("RF").cd(0);
        this.getAnalysisCanvas().getCanvas("RF").draw(this.getDataGroup().getItem(4).getH2F("hi_rf_pos_paddle"));
        this.getAnalysisCanvas().getCanvas("RF").cd(1);
        this.getAnalysisCanvas().getCanvas("RF").draw(this.getDataGroup().getItem(4).getH2F("hi_rf_neg_paddle"));
        this.getAnalysisCanvas().getCanvas("RF").cd(2);
        this.getAnalysisCanvas().getCanvas("RF").draw(this.getDataGroup().getItem(4).getH1F("hi_rf_pos"));
        this.getAnalysisCanvas().getCanvas("RF").draw(this.getDataGroup().getItem(4).getF1D("f1_rf_pos"),"same");
        this.getAnalysisCanvas().getCanvas("RF").cd(3);
        this.getAnalysisCanvas().getCanvas("RF").getPad(3).getAxisY().setRange(0, 0.3);
        this.getAnalysisCanvas().getCanvas("RF").draw(this.getDataGroup().getItem(4).getGraph("g_rf_pos_paddle"));
        this.getAnalysisCanvas().getCanvas("RF").cd(4);
        this.getAnalysisCanvas().getCanvas("RF").getPad(4).getAxisY().setRange(0, 0.3);
        this.getAnalysisCanvas().getCanvas("RF").draw(this.getDataGroup().getItem(4).getGraph("g_rf_neg_paddle"));
        this.getAnalysisCanvas().getCanvas("RF").cd(5);
        this.getAnalysisCanvas().getCanvas("RF").draw(this.getDataGroup().getItem(4).getH1F("hi_rf_neg"));
        this.getAnalysisCanvas().getCanvas("RF").draw(this.getDataGroup().getItem(4).getF1D("f1_rf_neg"),"same");
        this.getAnalysisCanvas().getCanvas("Residuals").cd(0);
        this.getAnalysisCanvas().getCanvas("Residuals").draw(this.getDataGroup().getItem(3).getH1F("hi_z_hit"));
        this.getAnalysisCanvas().getCanvas("Residuals").cd(1);
        this.getAnalysisCanvas().getCanvas("Residuals").draw(this.getDataGroup().getItem(3).getH2F("hi_rz_hit"));
        this.getAnalysisCanvas().getCanvas("Residuals").cd(2);
        this.getAnalysisCanvas().getCanvas("Residuals").draw(this.getDataGroup().getItem(3).getH2F("hi_z_hit_track"));
        this.getAnalysisCanvas().getCanvas("Residuals").cd(3);
        this.getAnalysisCanvas().getCanvas("Residuals").draw(this.getDataGroup().getItem(3).getH1F("hi_z_track"));
        this.getAnalysisCanvas().getCanvas("Residuals").cd(4);
        this.getAnalysisCanvas().getCanvas("Residuals").draw(this.getDataGroup().getItem(3).getH2F("hi_rz_track"));
        this.getAnalysisCanvas().getCanvas("Residuals").cd(5);
        this.getAnalysisCanvas().getCanvas("Residuals").getPad(5).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("Residuals").draw(this.getDataGroup().getItem(3).getH2F("hi_z_residual"));
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
        this.getAnalysisCanvas().getCanvas("Beta and Mass").draw(this.getDataGroup().getItem(5).getH1F("hi_mass_pos"));
        this.getAnalysisCanvas().getCanvas("Beta and Mass").cd(1);
        this.getAnalysisCanvas().getCanvas("Beta and Mass").getPad(1).getAxisY().setLog(true);
        this.getAnalysisCanvas().getCanvas("Beta and Mass").draw(this.getDataGroup().getItem(5).getH1F("hi_mass_neg"));
        this.getAnalysisCanvas().getCanvas("Beta and Mass").cd(2);
        this.getAnalysisCanvas().getCanvas("Beta and Mass").getPad(2).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("Beta and Mass").draw(this.getDataGroup().getItem(5).getH2F("hi_beta_pos"));
        this.getAnalysisCanvas().getCanvas("Beta and Mass").draw(this.getDataGroup().getItem(5).getF1D("cpion"),"same");
        this.getAnalysisCanvas().getCanvas("Beta and Mass").draw(this.getDataGroup().getItem(5).getF1D("ckaon"),"same");
        this.getAnalysisCanvas().getCanvas("Beta and Mass").draw(this.getDataGroup().getItem(5).getF1D("cprot"),"same");
        this.getAnalysisCanvas().getCanvas("Beta and Mass").cd(3);
        this.getAnalysisCanvas().getCanvas("Beta and Mass").getPad(3).getAxisZ().setLog(true);
        this.getAnalysisCanvas().getCanvas("Beta and Mass").draw(this.getDataGroup().getItem(5).getH2F("hi_beta_neg"));
        this.getAnalysisCanvas().getCanvas("Beta and Mass").draw(this.getDataGroup().getItem(5).getF1D("cpion"),"same");
        this.getAnalysisCanvas().getCanvas("Beta and Mass").draw(this.getDataGroup().getItem(5).getF1D("ckaon"),"same");

        this.getAnalysisCanvas().getCanvas("MIPs").update();    
        this.getAnalysisCanvas().getCanvas("ADC-TDC").update();    
        this.getAnalysisCanvas().getCanvas("RF").update();    
        this.getAnalysisCanvas().getCanvas("Residuals").update();    
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
        DataBank recTrackEB = null;
        DataBank recEvenEB = null;
        DataBank recRunRF  = null;
        DataBank recCndHits = null;
        DataBank recHBTTrack = null;
        DataBank cndADC = null;
        DataBank cndTDC = null;
        if(event.hasBank("RUN::config"))            recRun      = event.getBank("RUN::config");
        if(event.hasBank("REC::Particle"))          recBankEB   = event.getBank("REC::Particle");
        if(event.hasBank("REC::Scintillator"))      recDeteEB   = event.getBank("REC::Scintillator");
        if(event.hasBank("REC::Track"))             recTrackEB  = event.getBank("REC::Track");
        if(event.hasBank("REC::Event"))             recEvenEB   = event.getBank("REC::Event");
        if(event.hasBank("CND::hits"))              recCndHits = event.getBank("CND::hits");
        if(event.hasBank("CND::adc"))               cndADC     = event.getBank("CND::adc");
        if(event.hasBank("CND::tdc"))               cndTDC     = event.getBank("CND::tdc");
        if(event.hasBank("CVTRec::Tracks"))         recHBTTrack = event.getBank("CVTRec::Tracks");
        if(recRun == null) return;
        int ev  = recRun.getInt("event",0);
        int run = recRun.getInt("run",0);
        IndexedTable rfConfig = this.getCcdb().getConstants(run, "/calibration/eb/rf/config");
        if(run!=0 && this.rfPeriod!=rfConfig.getDoubleValue("clock", 1,1,1)) {
            this.rfPeriod = rfConfig.getDoubleValue("clock", 1,1,1);
            this.resetEventListener();
        }
//        System.out.println(ev); 
//        if(ev==134 || ev==370) {System.out.println(ev); recBankEB.show(); recDeteEB.show();}
//        if(ev==93364) recCndHits.show();
        if(recCndHits!=null) {
            int nrows = recCndHits.rows();
            for(int loop = 0; loop < nrows; loop++){
                int status    = recCndHits.getShort("status", loop);
                int sector    = recCndHits.getByte("sector", loop);
                int layer     = recCndHits.getByte("layer", loop);
                int comp      = recCndHits.getShort("component", loop);
                int paddle    = (layer-1)*48+(sector-1)*2+comp;
                int trk_id    = recCndHits.getShort("trkID", loop);
                double energy = recCndHits.getFloat("energy", loop);
                double time   = recCndHits.getFloat("time",loop);
           	double x      = recCndHits.getFloat("x", loop);
                double y      = recCndHits.getFloat("y", loop);
                double z      = recCndHits.getFloat("z", loop);
           	double tx     = recCndHits.getFloat("tx", loop);
                double ty     = recCndHits.getFloat("ty", loop);
                double tz     = recCndHits.getFloat("tz", loop);
                double path   = recCndHits.getFloat("pathlength", loop);
                double dx     = recCndHits.getFloat("tlength", loop);
        	int adcId1    = recCndHits.getShort("indexLadc", loop);
                int adcId2    = recCndHits.getShort("indexRadc", loop);
                int tdcId1    = recCndHits.getShort("indexLtdc", loop);
                int tdcId2    = recCndHits.getShort("indexRtdc", loop);
                double adc1   = 0;
                double adc2   = 0;
                double adct1  = 0;
                double adct2  = 0;
                double tdc1   = 0;
                double tdc2   = 0;
                if(cndADC!=null) {
                    adc1   = (double) cndADC.getInt("ADC",adcId1);
                    adc2   = (double) cndADC.getInt("ADC",adcId2);
                    adct1  = (double) cndADC.getFloat("time",adcId1);
                    adct2  = (double) cndADC.getFloat("time",adcId2);                    
                }
                if(cndTDC!=null) {
                    tdc1   = (double) cndTDC.getInt("TDC",tdcId1);
                    tdc2   = (double) cndTDC.getInt("TDC",tdcId2);                    
                }
                if(adc1>0 && adct1>0 && adc2>0 && adct2>0 ) {
                    this.getDataGroup().getItem(6).getH2F("hi_adctdc_left").fill(adct1,tdc1*tdcConv);
                    this.getDataGroup().getItem(6).getH2F("hi_adctdc_right").fill(adct2,tdc2*tdcConv);                
                    this.getDataGroup().getItem(6).getH1F("hi_adc_tdc_left").fill(tdc1*tdcConv-adct1);
                    this.getDataGroup().getItem(6).getH1F("hi_adc_tdc_right").fill(tdc2*tdcConv-adct2);                
                    this.getDataGroup().getItem(6).getH2F("hi_adc_tdc_paddle_left").fill(tdc1*tdcConv-adct1,paddle);
                    this.getDataGroup().getItem(6).getH2F("hi_adc_tdc_paddle_right").fill(tdc2*tdcConv-adct2,paddle); 
                }

                if(status>=0) {
                    this.getDataGroup().getItem(2).getH2F("hi_en_paddle").fill(paddle*1.,energy);
                    this.getDataGroup().getItem(2).getH2F("hi_time_paddle").fill(paddle*1.,(tdc1-tdc2)*tdcConv);
                    if(trk_id!=-1 && energy>0) {
                        int    q    = recHBTTrack.getInt("q",trk_id);
                        double p    = recHBTTrack.getFloat("p",trk_id);
                        double pt   = recHBTTrack.getFloat("pt",trk_id);
                        double phi0 = recHBTTrack.getFloat("phi0",trk_id);
                        Particle recParticle = new Particle(211,pt*Math.cos(phi0),pt*Math.sin(phi0),Math.sqrt(p*p-pt*pt),0,0,0);
                        double beta = recParticle.p()/Math.sqrt(recParticle.p()*recParticle.p()+recParticle.mass2());
                        this.getDataGroup().getItem(3).getH1F("hi_z_hit").fill(z);
                        this.getDataGroup().getItem(3).getH1F("hi_z_track").fill(tz);
                        this.getDataGroup().getItem(3).getH2F("hi_rz_hit").fill(z,Math.sqrt(x*x+y*y));
                        this.getDataGroup().getItem(3).getH2F("hi_rz_track").fill(tz,Math.sqrt(tx*tx+ty*ty));
                        this.getDataGroup().getItem(3).getH2F("hi_z_hit_track").fill(z,tz);
                        this.getDataGroup().getItem(3).getH2F("hi_z_residual").fill(z-tz,paddle*1.);
                        if(recEvenEB!=null && recBankEB!=null && recDeteEB!=null && recTrackEB!=null) {
                            double vt = -100;
                            int pindex = -1;
                            for(int i=0; i<recTrackEB.rows(); i++) {
                                if(recTrackEB.getShort("index",i)==trk_id && recTrackEB.getByte("detector",i)==DetectorType.CVT.getDetectorId()) {
                                    vt=recBankEB.getFloat("vt", recTrackEB.getShort("pindex", i));
                                    pindex=recTrackEB.getShort("pindex", i);
                                    break;
                                }                               
                            }
                            for(int i=0; i<recDeteEB.rows(); i++) {
                                if(recDeteEB.getShort("pindex",i)==pindex && recDeteEB.getByte("detector",i)==DetectorType.CND.getDetectorId()) {
                                    time=recDeteEB.getFloat("time", i);
                                    path=recDeteEB.getFloat("path", i);
                                    break;
                                }                               
                            }
                            if(vt>-100 && recBankEB.getInt("pid", 0)==11) {
                                double dt = time - path/(beta*29.97) - vt;
                                double betaTof = path/(time-vt)/29.97;
                                double mass2   = Math.pow(recParticle.p()/betaTof, 2)-recParticle.p()*recParticle.p();
                                if(q==1)  {
                                    this.getDataGroup().getItem(1).getH1F("hi_pos_en").fill(energy/dx);
                                    this.getDataGroup().getItem(1).getH2F("hi_pos_en_p").fill(recParticle.p(),energy/dx);
                                    this.getDataGroup().getItem(1).getH2F("hi_pos_en_paddle").fill(paddle*1.,energy/dx);
                                    if(recParticle.p()>0.4) {
                                        this.getDataGroup().getItem(4).getH2F("hi_rf_pos_paddle").fill(paddle*1.,dt);
                                        this.getDataGroup().getItem(4).getH1F("hi_rf_pos").fill(dt);
                                    }
                                    this.getDataGroup().getItem(5).getH2F("hi_beta_pos").fill(recParticle.p(),betaTof);
                                    this.getDataGroup().getItem(5).getH1F("hi_mass_pos").fill(mass2);
                                }
                                if(q==-1) {
                                    this.getDataGroup().getItem(1).getH1F("hi_neg_en").fill(energy/dx);
                                    this.getDataGroup().getItem(1).getH2F("hi_neg_en_p").fill(recParticle.p(),energy/dx);
                                    this.getDataGroup().getItem(1).getH2F("hi_neg_en_paddle").fill(paddle*1.,energy/dx);
                                    if(recParticle.p()>0.4) {
                                        this.getDataGroup().getItem(4).getH2F("hi_rf_neg_paddle").fill(paddle*1.,dt);
                                        this.getDataGroup().getItem(4).getH1F("hi_rf_neg").fill(dt);
                                    }
                                    this.getDataGroup().getItem(5).getH2F("hi_beta_neg").fill(recParticle.p(),betaTof);
                                    this.getDataGroup().getItem(5).getH1F("hi_mass_neg").fill(mass2);
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
//        System.out.println("Updating FTOF");
        this.analyze();
        //fitting negative tracks vertex
//        ParallelSliceFitter fitter = new ParallelSliceFitter(this.getDataGroup().getItem(5).getH2F("hi_rf_paddle_2"));
//        fitter.fitSlicesX();
//        GraphErrors meanY = fitter.getMeanSlices();
//        this.getDataGroup().getItem(5).getGraph("hi_rf_offsets_2").copy(meanY);       
    }
    
    @Override
    public void analyze() {
        this.fitSlices(this.getDataGroup().getItem(4).getH2F("hi_rf_pos_paddle"), this.getDataGroup().getItem(4).getGraph("g_rf_pos_paddle"));
        this.fitSlices(this.getDataGroup().getItem(4).getH2F("hi_rf_neg_paddle"), this.getDataGroup().getItem(4).getGraph("g_rf_neg_paddle"));
        this.initParGauss(this.getDataGroup().getItem(4).getH1F("hi_rf_pos"), this.getDataGroup().getItem(4).getF1D("f1_rf_pos"));
        this.initParGauss(this.getDataGroup().getItem(4).getH1F("hi_rf_neg"), this.getDataGroup().getItem(4).getF1D("f1_rf_neg"));
    }

    public void fitSlices(H2F histo, GraphErrors graph) {
        graph.reset();
        ArrayList<H1F> hslice = histo.getSlicesX();
        for(int i=0; i<hslice.size(); i++) {
            double  x = histo.getXAxis().getBinCenter(i);
            double ex = 0;
            double  y = hslice.get(i).getRMS();
            double ey = 0;
            F1D f1 = new F1D("f1","[amp]*gaus(x,[mean],[sigma])", -1.0, 1.0);
            this.initParGauss(hslice.get(i), f1);
            if(f1.getParameter(0)>50) graph.addPoint(x, f1.getParameter(2), ex, f1.parameter(2).error());
            else                      graph.addPoint(x, 0, ex, 0);
        }
    }
    
    public void initParGauss(H1F histo, F1D f1) {
        double mean  = histo.getDataX(histo.getMaximumBin());
        double amp   = histo.getBinContent(histo.getMaximumBin());
        double sigma = histo.getRMS();
        f1.setParameter(0, amp);
        f1.setParameter(1, mean);
        f1.setParameter(2, 0.1);
        f1.setRange(mean-1.6*sigma, mean+1.6*sigma);
        DataFitter.fit(f1, histo, "Q"); //No options uses error for sigma 
    }
}