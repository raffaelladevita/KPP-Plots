/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.analysis;

import org.clas.viewer.AnalysisMonitor;
import org.jlab.clas.physics.Particle;
import org.jlab.groot.data.H1F;
import org.jlab.groot.group.DataGroup;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;

/**
 *
 * @author devita
 */
public class TIMEmonitor extends AnalysisMonitor {
    

    public TIMEmonitor(String name) {
        super(name);
        this.setAnalysisTabNames("Timing Correlations");
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
        // electrons
        H1F hi_dc_1 = new H1F("hi_dc_1","hi_dc_1",200,-400,1000);
        hi_dc_1.setTitleX("Time (ns)");
        hi_dc_1.setTitleY("Counts");
        hi_dc_1.setTitle("DC region 1 ");
        H1F hi_dc_2 = new H1F("hi_dc_2","hi_dc_2",200,-400,1000);
        hi_dc_2.setTitleX("Time (ns)");
        hi_dc_2.setTitleY("Counts");
        hi_dc_2.setTitle("DC region 2 ");
        H1F hi_dc_3 = new H1F("hi_dc_3","hi_dc_3",200,-400,1000);
        hi_dc_3.setTitleX("Time (ns)");
        hi_dc_3.setTitleY("Counts");
        hi_dc_3.setTitle("DC region 3 ");
        H1F hi_ftof_1 = new H1F("hi_ftof_1","hi_ftof_1",200,400,1000);
        hi_ftof_1.setTitleX("Time (ns)");
        hi_ftof_1.setTitleY("Counts");
        hi_ftof_1.setTitle("FTOF panel 1A ");
        hi_ftof_1.setFillColor(23);
        H1F hi_ftof_2 = new H1F("hi_ftof_2","hi_ftof_2",200,400,1000);
        hi_ftof_2.setTitleX("Time (ns)");
        hi_ftof_2.setTitleY("Counts");
        hi_ftof_2.setTitle("FTOF panel 1B ");
        hi_ftof_2.setFillColor(23);
        H1F hi_ftof_3 = new H1F("hi_ftof_3","hi_ftof_3",200,400,1000);
        hi_ftof_3.setTitleX("Time (ns)");
        hi_ftof_3.setTitleY("Counts");
        hi_ftof_3.setTitle("FTOF panel 2 ");
        hi_ftof_3.setFillColor(23);
        H1F hi_ecal_1 = new H1F("hi_ecal_1","hi_ecal_1",200,400,1000);
        hi_ecal_1.setTitleX("Time (ns)");
        hi_ecal_1.setTitleY("Counts");
        hi_ecal_1.setTitle("PCAL ");
        H1F hi_ecal_2 = new H1F("hi_ecal_2","hi_ecal_2",200,400,1000);
        hi_ecal_2.setTitleX("Time (ns)");
        hi_ecal_2.setTitleY("Counts");
        hi_ecal_2.setTitle("EC inner ");
        H1F hi_ecal_3 = new H1F("hi_ecal_3","hi_ecal_3",200,400,1000);
        hi_ecal_3.setTitleX("Time (ns)");
        hi_ecal_3.setTitleY("Counts");
        hi_ecal_3.setTitle("EC outer ");
        DataGroup dg_time = new DataGroup(3,3);
        dg_time.addDataSet(hi_dc_1, 0);
        dg_time.addDataSet(hi_dc_2, 3);
        dg_time.addDataSet(hi_dc_3, 6);
        dg_time.addDataSet(hi_ftof_1, 1);
        dg_time.addDataSet(hi_ftof_2, 4);
        dg_time.addDataSet(hi_ftof_3, 7);
        dg_time.addDataSet(hi_ecal_1, 2);
        dg_time.addDataSet(hi_ecal_2, 5);
        dg_time.addDataSet(hi_ecal_3, 8);
        this.getDataGroup().add(dg_time, 1);       
    }
    
    @Override
    public void plotHistos() {
        this.getAnalysisCanvas().getCanvas("Timing Correlations").divide(1, 3);
        this.getAnalysisCanvas().getCanvas("Timing Correlations").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Timing Correlations").setGridY(false);
        
        for(int i=0; i<3; i++) {
            String hname = "hi_ftof_" + (i+1);
            this.getAnalysisCanvas().getCanvas("Timing Correlations").cd(i);
            this.getAnalysisCanvas().getCanvas("Timing Correlations").getPad(i).setTitleFontSize(24);
            this.getAnalysisCanvas().getCanvas("Timing Correlations").draw(this.getDataGroup().getItem(1).getH1F(hname));
        }

        this.getAnalysisCanvas().getCanvas("Timing Correlations").update();
    }
        
    @Override
    public void processEvent(DataEvent event) {
        // process event info and save into data group        
        Particle partRecEB = null;

        // event builder
        DataBank ftof = event.getBank("FTOF::tdc");
        DataBank ecal = event.getBank("ECAL::tdc");
        DataBank dc   = event.getBank("DC::tdc");
        if(ftof!=null) {
            int nrows = ftof.rows();
            int rows = ftof.rows();
            for(int i = 0; i < rows; i++){
                int    sector = ftof.getByte("sector",i);
                int     layer = ftof.getByte("layer",i);
                int    paddle = ftof.getShort("component",i);
                int       TDC = ftof.getInt("TDC",i);
                int     order = ftof.getByte("order",i); // order specifies left-right for ADC
                this.getDataGroup().getItem(1).getH1F("hi_ftof_" + layer).fill(TDC*0.0248);
            }
        }
        if(ecal!=null) {
            int nrows = ecal.rows();
            int rows = ecal.rows();
            for(int i = 0; i < rows; i++){
                int    sector = ecal.getByte("sector",i);
                int     layer = ecal.getByte("layer",i);
                int    paddle = ecal.getShort("component",i);
                int       TDC = ecal.getInt("TDC",i);
                int     order = ecal.getByte("order",i); // order specifies left-right for ADC
                int superlayer = ((int) (layer-1)/3)+1;
                this.getDataGroup().getItem(1).getH1F("hi_ecal_" + superlayer).fill(TDC*0.0248);
            }
        }
        if(ecal!=null) {
            int nrows = dc.rows();
            int rows = dc.rows();
            for(int i = 0; i < rows; i++){
                int    sector = dc.getByte("sector",i);
                int     layer = dc.getByte("layer",i);
                int      wire = dc.getShort("component",i);
                int       TDC = dc.getInt("TDC",i);
                int     order = dc.getByte("order",i); // order specifies left-right for ADC
                int    region = ((int) (layer-1)/12)+1;
                this.getDataGroup().getItem(1).getH1F("hi_dc_" + region).fill(TDC*0.1);
            }
        }
    }
    
    @Override
    public void timerUpdate() {
//        System.out.println("Updating TIME");
    }
}
