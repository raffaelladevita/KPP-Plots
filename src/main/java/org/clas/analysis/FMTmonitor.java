/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.analysis;

import java.util.ArrayList;
import java.util.List;
import org.clas.viewer.AnalysisMonitor;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.group.DataGroup;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;

/**
 *
 * @author devita
 */
public class FMTmonitor extends AnalysisMonitor {

    double rfPeriod = 4.008;
    public static final int FVT_Nlayers = 6;
    public static int FVT_Nstrips = 1024;// Number of strips
    public static int FVT_Halfstrips = 320; // In the middle of the FMT, 320 strips are split in two. 
    public static double FVT_Pitch = 0.0525; //strip width
    public static double FVT_Beamhole = 4.2575;//Radius of the hole in the center for the beam.

    double[] FVT_stripsXlocref;
    double[] FVT_stripsYlocref;
    double[] FVT_stripslength; //Give the strip length
    double[][][] FVT_stripsX; //Give the  end-points x-coordinates of the strip segment rotated in the correct frame for the layer
    double[][][] FVT_stripsY; //Give the  end-points y-coordinates of the strip segment

    double[] FVT_Zlayer = {294.997, 306.897, 318.797, 332.697, 344.597, 356.497}; //Give z-coordinate of the layer
    double[] FVT_Alpha = {19, 79, 139, -161, -101, -41}; //Give the rotation angle to apply
    double[] MY_Alpha = {170, 50, -70, 170, 50, -70}; //Give the rotation angle to apply

    public FMTmonitor(String name) {
        super(name);
        this.setAnalysisTabNames("FMT hits", "FMT phi", "Clusters", "Cluster phi", "Cluster phi vs strip", "Cluster strip vs strip", "Cluster energy");
        this.init(false);
        this.LoadGeometry();
    }

    @Override
    public void createHistos() {
        // initialize canvas and create histograms
        this.setNumberOfEvents(0);
        H1F summary = new H1F("summary", "summary", 6, 1, 7);
        summary.setTitleX("sector");
        summary.setTitleY("DC hits");
        summary.setFillColor(33);
        DataGroup sum = new DataGroup(1, 1);
        sum.addDataSet(summary, 0);
        this.setAnalysisSummary(sum);

        // Hits
        DataGroup dc_hits = new DataGroup(3, 4);
        for (int layer = 1; layer <= FVT_Nlayers; layer++) {
            H1F hi_hit_res = new H1F("hi_hit_res_l" + layer, "Residual (cm) - Layer " + layer, "Counts", 200, -30, 30);
            hi_hit_res.setFillColor(4);
            H1F hi_hit_phi = new H1F("hi_hit_phi_l" + layer, "Strip angle - Layer " + layer, "Counts", 200, -180, 180);
            hi_hit_phi.setFillColor(4);
            H1F hi_hit_phi_cut = new H1F("hi_hit_phi_cut_l" + layer, "Strip angle - Layer " + layer, "Counts", 200, -180, 180);
            hi_hit_phi_cut.setFillColor(3);
            dc_hits.addDataSet(hi_hit_res, layer - 1);
            dc_hits.addDataSet(hi_hit_phi, layer - 1 + FVT_Nlayers);
            dc_hits.addDataSet(hi_hit_phi_cut, layer - 1 + FVT_Nlayers);
        }
        this.getDataGroup().add(dc_hits, 0);
        // Clusters
        DataGroup dc_clusters = new DataGroup(3, 10);
        for (int layer = 1; layer <= FVT_Nlayers; layer++) {
            H1F hi_cluster_res = new H1F("hi_cluster_res_l" + layer, "Residual (cm) - Layer " + layer, "Counts", 200, -30, 30);
            hi_cluster_res.setFillColor(4);
            H1F hi_cluster_phi = new H1F("hi_cluster_phi_l" + layer, "Strip angle - Layer " + layer, "Counts", 200, -180, 180);
            hi_cluster_phi.setFillColor(4);
            H1F hi_cluster_phi_cut = new H1F("hi_cluster_phi_cut_l" + layer, "Strip angle - Layer " + layer, "Counts", 200, -180, 180);
            hi_cluster_phi_cut.setFillColor(3);
            H1F hi_cluster_energy_cut = new H1F("hi_cluster_energy_cut_l" + layer, "Strip angle - Layer " + layer, "Counts", 200, -180, 180);
            hi_cluster_energy_cut.setFillColor(2);
            H2F hi_cluster_phi_energy = new H2F("hi_cluster_phi_energy_l" + layer, 100, -180, 180, 100, 0, 500);
            hi_cluster_phi_energy.setTitleX("Strip angle - Layer " + layer);
            hi_cluster_phi_energy.setTitleY("Energy");
            H2F hi_cluster_phi_strip = new H2F("hi_cluster_phi_strip_l" + layer, 100, -180, 180, 100, 0, 1000);
            hi_cluster_phi_strip.setTitleX("Strip angle - Layer " + layer);
            hi_cluster_phi_strip.setTitleY("Strip - Layer " + layer);
            H2F hi_cluster_strip_strip = new H2F("hi_cluster_strip_strip_l" + layer, 100, 0, 1000, 100, 0, 1000);
            hi_cluster_strip_strip.setTitleX("Strip - Layer " + "3");
            hi_cluster_strip_strip.setTitleY("Strip - Layer " + layer);
            dc_clusters.addDataSet(hi_cluster_res, layer - 1);
            dc_clusters.addDataSet(hi_cluster_phi, layer - 1 + FVT_Nlayers);
            dc_clusters.addDataSet(hi_cluster_phi_cut, layer - 1 + FVT_Nlayers);
            dc_clusters.addDataSet(hi_cluster_energy_cut, layer - 1 + FVT_Nlayers);
            dc_clusters.addDataSet(hi_cluster_phi_strip, layer - 1 + 2*FVT_Nlayers);
            dc_clusters.addDataSet(hi_cluster_strip_strip, layer - 1 + 3*FVT_Nlayers);
            dc_clusters.addDataSet(hi_cluster_phi_energy, layer - 1 + 4*FVT_Nlayers);
        }
        this.getDataGroup().add(dc_clusters, 1);

    }

    @Override
    public void plotHistos() {

        this.getAnalysisCanvas().getCanvas("FMT hits").divide(3, 2);
        this.getAnalysisCanvas().getCanvas("FMT hits").setGridX(false);
        this.getAnalysisCanvas().getCanvas("FMT hits").setGridY(false);
        this.getAnalysisCanvas().getCanvas("FMT phi").divide(3, 2);
        this.getAnalysisCanvas().getCanvas("FMT phi").setGridX(false);
        this.getAnalysisCanvas().getCanvas("FMT phi").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Clusters").divide(3, 2);
        this.getAnalysisCanvas().getCanvas("Clusters").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Clusters").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Cluster phi").divide(3, 2);
        this.getAnalysisCanvas().getCanvas("Cluster phi").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Cluster phi").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Cluster phi vs strip").divide(3, 2);
        this.getAnalysisCanvas().getCanvas("Cluster phi vs strip").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Cluster phi vs strip").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Cluster strip vs strip").divide(3, 2);
        this.getAnalysisCanvas().getCanvas("Cluster strip vs strip").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Cluster strip vs strip").setGridY(false);
        this.getAnalysisCanvas().getCanvas("Cluster energy").divide(3, 2);
        this.getAnalysisCanvas().getCanvas("Cluster energy").setGridX(false);
        this.getAnalysisCanvas().getCanvas("Cluster energy").setGridY(false);

        for (int layer = 1; layer <= FVT_Nlayers; layer++) {
            this.getAnalysisCanvas().getCanvas("FMT hits").cd(layer - 1);
//        this.getAnalysisCanvas().getCanvas("Beta").getPad(0).getAxisZ().setLog(true);
            this.getAnalysisCanvas().getCanvas("FMT hits").draw(this.getDataGroup().getItem(0).getH1F("hi_hit_res_l" + layer));
            this.getAnalysisCanvas().getCanvas("FMT phi").cd(layer - 1);
            this.getAnalysisCanvas().getCanvas("FMT phi").draw(this.getDataGroup().getItem(0).getH1F("hi_hit_phi_l" + layer));
            this.getAnalysisCanvas().getCanvas("FMT phi").draw(this.getDataGroup().getItem(0).getH1F("hi_hit_phi_cut_l" + layer), "same");
        }
        for (int layer = 1; layer <= FVT_Nlayers; layer++) {
            this.getAnalysisCanvas().getCanvas("Clusters").cd(layer - 1);
//        this.getAnalysisCanvas().getCanvas("Beta").getPad(0).getAxisZ().setLog(true);
            this.getAnalysisCanvas().getCanvas("Clusters").draw(this.getDataGroup().getItem(1).getH1F("hi_cluster_res_l" + layer));
            this.getAnalysisCanvas().getCanvas("Cluster phi").cd(layer - 1);
            this.getAnalysisCanvas().getCanvas("Cluster phi").draw(this.getDataGroup().getItem(1).getH1F("hi_cluster_phi_l" + layer));
            this.getAnalysisCanvas().getCanvas("Cluster phi").draw(this.getDataGroup().getItem(1).getH1F("hi_cluster_energy_cut_l" + layer), "same");
            this.getAnalysisCanvas().getCanvas("Cluster phi").draw(this.getDataGroup().getItem(1).getH1F("hi_cluster_phi_cut_l" + layer), "same");
            this.getAnalysisCanvas().getCanvas("Cluster phi vs strip").cd(layer - 1);
//            this.getAnalysisCanvas().getCanvas("Cluster phi vs strip").getPad(layer - 1).getAxisZ().setLog(true);
            this.getAnalysisCanvas().getCanvas("Cluster phi vs strip").draw(this.getDataGroup().getItem(1).getH2F("hi_cluster_phi_strip_l" + layer));
            this.getAnalysisCanvas().getCanvas("Cluster strip vs strip").cd(layer - 1);
            this.getAnalysisCanvas().getCanvas("Cluster strip vs strip").getPad(layer - 1).getAxisZ().setLog(true);
            this.getAnalysisCanvas().getCanvas("Cluster strip vs strip").draw(this.getDataGroup().getItem(1).getH2F("hi_cluster_strip_strip_l" + layer));
            this.getAnalysisCanvas().getCanvas("Cluster energy").cd(layer - 1);
            this.getAnalysisCanvas().getCanvas("Cluster energy").draw(this.getDataGroup().getItem(1).getH2F("hi_cluster_phi_energy_l" + layer));
        }
    }

    @Override
    public void processEvent(DataEvent event) {
        // process event info and save into data group        
        // event builder
        DataBank recRun = null;
        DataBank recBankEB = null;
        DataBank recEvenEB = null;
        DataBank recTrajEB = null;
        DataBank fmtHits = null;
        DataBank fmtClusters = null;
        if (event.hasBank("RUN::config")) {
            recRun = event.getBank("RUN::config");
        }
        if (event.hasBank("REC::Particle")) {
            recBankEB = event.getBank("REC::Particle");
        }
        if (event.hasBank("REC::Event")) {
            recEvenEB = event.getBank("REC::Event");
        }
        if (event.hasBank("REC::Traj")) {
            recTrajEB = event.getBank("REC::Traj");
        }
        if (event.hasBank("FMTRec::Hits")) {
            fmtHits = event.getBank("FMTRec::Hits");
        }
        if (event.hasBank("FMTRec::Clusters")) {
            fmtClusters = event.getBank("FMTRec::Clusters");
        }
        int ev = 0;
        if(recRun != null) ev=recRun.getInt("event",0);
//        System.out.println(ev); 

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
        double startTime = -1000;
        double rfTime = -1000;
        if (recEvenEB != null) {
            startTime = recEvenEB.getFloat("STTime", 0);
            rfTime = recEvenEB.getFloat("RFTime", 0);
        }
        // get trigger particle
        int trigger = 0;
        if (recBankEB != null) {
            trigger = recBankEB.getInt("pid", 0);
        }

        if (fmtHits != null && fmtClusters != null && recTrajEB != null) {
//            System.out.println("Updating AAA");
            int[] matchedLayers = {0, 0, 0, 0, 0, 0};
            List<Integer> matchedStrips = new ArrayList<Integer>();
            List<Integer> matchedClusters = new ArrayList<Integer>();
            for (int loop = 0; loop < recTrajEB.rows(); loop++) {
                int detId = recTrajEB.getShort("detId", loop);
                if (detId >= 1 && detId <= FVT_Nlayers) {
                    double x = recTrajEB.getFloat("x", loop);
                    double y = recTrajEB.getFloat("y", loop);
                    double z = recTrajEB.getFloat("z", loop);
                    double phiRef = 0;
                    if (event.hasBank("MC::Particle")) {
                        phiRef = FVT_Alpha[detId - 1];
                    } else {
                        phiRef = MY_Alpha[detId - 1];
                    }
                    for (int i = 0; i < fmtHits.rows(); i++) {
                        if (detId == fmtHits.getByte("layer", i)) {
                            int strip = fmtHits.getInt("strip", i);
                            double xLoc = x * Math.cos(Math.toRadians(phiRef)) + y * Math.sin(Math.toRadians(phiRef));
                            double yLoc = y * Math.cos(Math.toRadians(phiRef)) - x * Math.sin(Math.toRadians(phiRef));
                            double phi[] = this.getPhi(x, y, FVT_stripsYlocref[strip]);
//                            if (Math.abs(yLoc - FVT_stripsYlocref[strip]) < 5) {
//                                matchedLayers[detId - 1]++;
//                            }
                            if(detId==3) matchedStrips.add(strip);
                            this.getDataGroup().getItem(0).getH1F("hi_hit_res_l" + detId).fill(yLoc - FVT_stripsYlocref[strip]);
                            if (phi[0] > -9999 && strip<=512) {
                                this.getDataGroup().getItem(0).getH1F("hi_hit_phi_l" + detId).fill(phi[0]);
                            }
                            if (phi[1] > -9999 && strip>512) {
                                this.getDataGroup().getItem(0).getH1F("hi_hit_phi_l" + detId).fill(phi[1]);
                            }
                        }
                    }
                    for (int i = 0; i < fmtClusters.rows(); i++) {
                        if (detId == fmtClusters.getByte("layer", i)) {
                            int strip = fmtClusters.getInt("seedStrip", i);
                            double energy = fmtClusters.getFloat("ETot", i);
                            double xLoc = x * Math.cos(Math.toRadians(phiRef)) + y * Math.sin(Math.toRadians(phiRef));
                            double yLoc = y * Math.cos(Math.toRadians(phiRef)) - x * Math.sin(Math.toRadians(phiRef));
                            double phi[] = this.getPhi(x, y, FVT_stripsYlocref[strip]);
//                            if (Math.abs(yLoc - FVT_stripsYlocref[strip]) < 5) {
//                                matchedLayers[detId - 1]++;
//                            }
                            if(detId==3) matchedClusters.add(strip);
                            this.getDataGroup().getItem(1).getH1F("hi_cluster_res_l" + detId).fill(yLoc - FVT_stripsYlocref[strip]);
                            if (phi[0] > -9999 && strip<=512) {
                                this.getDataGroup().getItem(1).getH1F("hi_cluster_phi_l" + detId).fill(phi[0]);
                                this.getDataGroup().getItem(1).getH2F("hi_cluster_phi_energy_l" + detId).fill(phi[0],energy);
                                if(energy>100) this.getDataGroup().getItem(1).getH1F("hi_cluster_energy_cut_l" + detId).fill(phi[0]);
                                this.getDataGroup().getItem(1).getH2F("hi_cluster_phi_strip_l" + detId).fill(phi[0],strip);
                            }
                            if (phi[1] > -9999 && strip>512) {
                                this.getDataGroup().getItem(1).getH1F("hi_cluster_phi_l" + detId).fill(phi[1]);
                                this.getDataGroup().getItem(1).getH2F("hi_cluster_phi_energy_l" + detId).fill(phi[1],energy);
                                if(energy>100) this.getDataGroup().getItem(1).getH1F("hi_cluster_energy_cut_l" + detId).fill(phi[1]);
                                this.getDataGroup().getItem(1).getH2F("hi_cluster_phi_strip_l" + detId).fill(phi[1],strip);
                            }
                        }
                    }
                }
            }
            for (int loop = 0; loop < recTrajEB.rows(); loop++) {
                int detId = recTrajEB.getShort("detId", loop);
                if (detId >= 1 && detId <= FVT_Nlayers) {
                    double x = recTrajEB.getFloat("x", loop);
                    double y = recTrajEB.getFloat("y", loop);
                    double z = recTrajEB.getFloat("z", loop);
//                    if (matchedLayers[0] > 0 && matchedLayers[1] > 0 && matchedLayers[3] > 0 && matchedLayers[4] > 0) {
                    for (int i = 0; i < fmtHits.rows(); i++) {
                        if (detId == fmtHits.getByte("layer", i)) {
                            int strip = fmtHits.getInt("strip", i);
                            boolean match = false;
                            for(int strip3 : matchedStrips) if(detId==6 && Math.abs(strip-strip3)<30) match=true;
                            if(match) {
                                double phi[] = this.getPhi(x, y, FVT_stripsYlocref[strip]);
                                if (phi[0] > -9999 && strip<=512) {
                                    this.getDataGroup().getItem(0).getH1F("hi_hit_phi_cut_l" + detId).fill(phi[0]);
                                }
                                if (phi[1] > -9999 && strip>512) {
                                    this.getDataGroup().getItem(0).getH1F("hi_hit_phi_cut_l" + detId).fill(phi[1]);
                                }
                            }
                        }
                    }
                    for (int i = 0; i < fmtClusters.rows(); i++) {
                        if (detId == fmtClusters.getByte("layer", i)) {
                            int strip = fmtClusters.getInt("seedStrip", i);
                            double energy = fmtClusters.getFloat("ETot", i);
                            boolean match = false;
                            for(int strip3 : matchedClusters) {
                                if(detId==6 && Math.abs(strip-strip3)<30) match=true;
                                this.getDataGroup().getItem(1).getH2F("hi_cluster_strip_strip_l" + detId).fill(strip3,strip);
                            }
                            if(match) {
                                double phi[] = this.getPhi(x, y, FVT_stripsYlocref[strip]);
                                if (phi[0] > -9999 && strip<=512 && energy>100) {
                                    this.getDataGroup().getItem(1).getH1F("hi_cluster_phi_cut_l" + detId).fill(phi[0]);
                                }
                                if (phi[1] > -9999 && strip>512 && energy>100) {
                                    this.getDataGroup().getItem(1).getH1F("hi_cluster_phi_cut_l" + detId).fill(phi[1]);
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
//        ParallelSliceFitter fitter = new ParallelSliceFitter(this.getDataGroup().getItem(4).getH2F("hi_rf_paddle_2"));
//        fitter.fitSlicesX();
//        GraphErrors meanY = fitter.getMeanSlices();
//        this.getDataGroup().getItem(4).getGraph("hi_rf_offsets_2").copy(meanY);       
    }

    @Override
    public void analyze() {
        // fit hodoscope charge
    }

    public void LoadGeometry() {

        int FVT_Sidestrips = (FVT_Nstrips - 2 * FVT_Halfstrips) / 2;
        double FVT_YCentral = (double) FVT_Halfstrips * FVT_Pitch / 2.;
        double FVT_Rmax = FVT_Pitch * (FVT_Halfstrips + 2 * FVT_Sidestrips) / 2.;
        double FVT_SigmaS = FVT_Pitch / Math.sqrt(12);
        double[][] FVT_stripsXloc = new double[FVT_Nstrips][2];
        double[][] FVT_stripsYloc = new double[FVT_Nstrips][2];
        FVT_stripsXlocref = new double[FVT_Nstrips];
        FVT_stripsYlocref = new double[FVT_Nstrips];
        FVT_stripsX = new double[FVT_Nlayers][FVT_Nstrips][2];
        FVT_stripsY = new double[FVT_Nlayers][FVT_Nstrips][2];
        FVT_stripslength = new double[FVT_Nstrips];

        for (int i = 0; i < FVT_Nstrips; i++) {
            //Give the Y of the middle of the strip
            if (i<512){
                    FVT_stripsYloc[i][0]=-FVT_Rmax+(511-i+0.5)*FVT_Pitch;
                    FVT_stripsYloc[i][1]=-FVT_Rmax+(511-i+0.5)*FVT_Pitch;
            } else {
                    FVT_stripsYloc[i][0]=FVT_Rmax-(1023-i+0.5)*FVT_Pitch;
                    FVT_stripsYloc[i][1]=FVT_Rmax-(1023-i+0.5)*FVT_Pitch;
            }
            FVT_stripsYlocref[i] = FVT_stripsYloc[i][0];

            int localRegion = getLocalRegion(i);
            switch (localRegion) {
                case 2:
                case 4:
                    FVT_stripslength[i] = 2 * FVT_Rmax * Math.sin(Math.acos(Math.abs(FVT_stripsYloc[i][0]) / FVT_Rmax));
                    FVT_stripsXloc[i][0] = -FVT_stripslength[i] / 2.;
                    FVT_stripsXloc[i][1] = FVT_stripslength[i] / 2.;
                    FVT_stripsXlocref[i] = 0;
                    break;
                case 1:
                    FVT_stripslength[i] = FVT_Rmax * Math.sin(Math.acos(Math.abs(FVT_stripsYloc[i][0]) / FVT_Rmax));
                    FVT_stripsXloc[i][1] = 0;
                    FVT_stripsXloc[i][0] = -FVT_stripslength[i];
                    FVT_stripsXlocref[i] = -FVT_stripslength[i] / 2;
                    if (Math.abs(FVT_stripsYloc[i][0]) / FVT_Beamhole < 1) {
                        FVT_stripslength[i] = FVT_Rmax * Math.sin(Math.acos(Math.abs(FVT_stripsYloc[i][0]) / FVT_Rmax)) - FVT_Beamhole * Math.sin(Math.acos(Math.abs(FVT_stripsYloc[i][0]) / FVT_Beamhole));
                        FVT_stripsXloc[i][1] = -FVT_Beamhole * Math.sin(Math.acos(Math.abs(FVT_stripsYloc[i][0]) / FVT_Beamhole));
                        FVT_stripsXloc[i][0] = -FVT_stripslength[i];
                        FVT_stripsXlocref[i] = -FVT_stripslength[i] / 2 - FVT_Beamhole * Math.sin(Math.acos(Math.abs(FVT_stripsYloc[i][0]) / FVT_Beamhole));
                    }
                    break;
                case 3:
                    FVT_stripslength[i] = FVT_Rmax * Math.sin(Math.acos(Math.abs(FVT_stripsYloc[i][0]) / FVT_Rmax));
                    FVT_stripsXloc[i][0] = 0;
                    FVT_stripsXloc[i][1] = FVT_stripslength[i];
                    FVT_stripsXlocref[i] = FVT_stripslength[i] / 2;
                    if (Math.abs(FVT_stripsYloc[i][0]) / FVT_Beamhole < 1) {
                        FVT_stripslength[i] = FVT_Rmax * Math.sin(Math.acos(Math.abs(FVT_stripsYloc[i][0]) / FVT_Rmax)) - FVT_Beamhole * Math.sin(Math.acos(Math.abs(FVT_stripsYloc[i][0]) / FVT_Beamhole));
                        FVT_stripsXloc[i][0] = FVT_Beamhole * Math.sin(Math.acos(Math.abs(FVT_stripsYloc[i][0]) / FVT_Beamhole));
                        FVT_stripsXloc[i][1] = FVT_stripslength[i];
                        FVT_stripsXlocref[i] = FVT_stripslength[i] / 2 + FVT_Beamhole * Math.sin(Math.acos(Math.abs(FVT_stripsYloc[i][0]) / FVT_Beamhole));
                    }
                    break;
            }
            for (int j = 0; j < FVT_Nlayers; j++) { //X sign flip
                FVT_stripsX[j][i][0] = -(FVT_stripsXloc[i][0] * Math.cos(FVT_Alpha[j]) + FVT_stripsYloc[i][0] * Math.sin(FVT_Alpha[j]));
                FVT_stripsY[j][i][0] = -FVT_stripsXloc[i][0] * Math.sin(FVT_Alpha[j]) + FVT_stripsYloc[i][0] * Math.cos(FVT_Alpha[j]);
                FVT_stripsX[j][i][1] = -(FVT_stripsXloc[i][1] * Math.cos(FVT_Alpha[j]) + FVT_stripsYloc[i][1] * Math.sin(FVT_Alpha[j]));
                FVT_stripsY[j][i][1] = -FVT_stripsXloc[i][1] * Math.sin(FVT_Alpha[j]) + FVT_stripsYloc[i][1] * Math.cos(FVT_Alpha[j]);
            }

            //System.out.println(Constants.getLocalRegion(i)+" strip-1 = "+i+" x' "+FVT_stripsXloc[i][1]+" y' "+FVT_stripsYloc[i][1]+" length "+FVT_stripslength[i]+" FVT_Beamhole "+FVT_Beamhole);
        }

    }

    private int getLocalRegion(int i) {
		// To represent the geometry we divide the barrel micromega disk into 3 regions according to the strip numbering system.
        // Here i = strip_number -1;
        // Region 1 is the region in the negative x part of inner region: the strips range is from   1 to 320  (   0 <= i < 320)
        // Region 2 is the region in the negative y part of outer region: the strips range is from 321 to 512  ( 320 <= i < 512)
        // Region 3 is the region in the positive x part of inner region: the strips range is from 513 to 832  ( 512 <= i < 832)
        // Region 4 is the region in the positive y part of outer region: the strips range is from 833 to 1024 ( 832 <= i < 1024)

        int region = 0;
        if (i >= 0 && i < 320) {
            region = 1;
        }
        if (i >= 320 && i < 512) {
            region = 2;
        }
        if (i >= 512 && i < 832) {
            region = 3;
        }
        if (i >= 832 && i < 1024) {
            region = 4;
        }

        return region;
    }

    private double[] getPhi(double x, double y, double yloc) {
        double t1 = -9999;
        double t2 = -9999;
        if (y == -yloc) {
            t1 = y / x;
            t2 = y / x;
        } else if (x * x > (yloc * yloc - y * y)) {
            t1 = (-x + Math.sqrt(x * x - (yloc * yloc - y * y))) / (yloc + y);
            t2 = (-x - Math.sqrt(x * x - (yloc * yloc - y * y))) / (yloc + y);
//            if(Math.abs(yloc-y*(1-t1*t1)/(1+t1*t1))-x*2*t1/(1+t1*t1)<0.1 || true) t=t1;
//            else t=t2;
        }
        double phi[] = new double[2];
        phi[0]=2 * Math.toDegrees(Math.atan(t1));
        phi[1]=2 * Math.toDegrees(Math.atan(t2));
        if(t1==-9999) phi[0]=-9999;
        if(t2==-9999) phi[1]=-9999;
        return phi;
    }
}
