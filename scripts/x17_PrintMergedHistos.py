import ROOT
import sys
import json

from ROOT import kBlack
from ROOT import kBlue
from ROOT import kMagenta
from ROOT import kPink
from ROOT import kRed
from ROOT import kGreen
from ROOT import kCyan
from ROOT import kSpring
from ROOT import gPad
from ROOT import gStyle
from pathlib import Path

ROOT.gErrorIgnoreLevel = ROOT.kWarning

def extract_histograms(directory, path=""):
    """Recursively traverses a ROOT directory to find all histograms."""
    histograms = {}
    
    # Loop over all keys inside the current directory
    for key in directory.GetListOfKeys():
        name = key.GetName()
        classname = key.GetClassName()
        full_path = f"{path}/{name}" if path else name
        
        # Retrieve the object from the file
        obj = key.ReadObj()
        
        # Check if the object is a subdirectory
        if isinstance(obj, ROOT.TDirectoryFile):
            # Recurse into the subdirectory
            histograms.update(extract_histograms(obj, full_path))
            
        # Check if the object is a histogram (inherits from TH1)
        elif obj.InheritsFrom("TH1"):
            # Detach histogram from the directory to keep it in memory after file closes
            obj.SetDirectory(0)  
            histograms[full_path] = obj
            print(f"Found histogram [{classname}]: {full_path}")
            
    return histograms

def printXPdf(X_histos, X_pdfName, lC, fo):
    #Separate out the histograms for the X candidates
    xy_histos = [w for w in X_histos if "HC_XY_" in w.GetName()]
    xy_ProjX = [w.ProjectionX() for w in xy_histos]
    xy_ProjY = [w.ProjectionY() for w in xy_histos]

    ETheta_histos = [w for w in X_histos if "E_theta_" in w.GetName()]
    sumPt_histos = [w for w in X_histos if "Sum_pt_" in w.GetName()]
    sumPxVPy_histos = [w for w in X_histos if "_pxVpy_" in w.GetName()]
    PxvPy_ProjX = [w.ProjectionX() for w in sumPxVPy_histos]
    PxvPy_ProjY = [w.ProjectionY() for w in sumPxVPy_histos]

    diffPhi_histos = [w for w in X_histos if "diffPhi_" in w.GetName()]
    time_histos = [w for w in X_histos if "_timing_" in w.GetName()]
    sumE_histos = [w for w in X_histos if "_sumE_" in w.GetName()]
    minE_histos = [w for w in X_histos if "_minE_" in w.GetName()]
    medE_histos = [w for w in X_histos if "_medE_" in w.GetName()]
    maxE_histos = [w for w in X_histos if "_maxE_" in w.GetName()]
    
    #Separate out the vertex Z of histograms for the X candidates
    vZ_histos = [w for w in X_histos if "_vZ_" in w.GetName()]
    vZ_All_histos = [w for w in vZ_histos if "_All_" in w.GetName()]
    vZ_MinE_histos = [w for w in vZ_histos if "_Min_E_" in w.GetName()]
    vZ_MedE_histos = [w for w in vZ_histos if "_Med_E_" in w.GetName()]
    vZ_MaxE_histos = [w for w in vZ_histos if "_Max_E_" in w.GetName()]

    mass_histos = [w for w in X_histos if "_invM_" in w.GetName()]
    lowBin = [mass_histos[0].GetXaxis().FindBin(16.75), mass_histos[0].GetXaxis().FindBin(16.25), mass_histos[0].GetXaxis().FindBin(15.75), mass_histos[0].GetXaxis().FindBin(14.75)]
    highBin = [mass_histos[0].GetXaxis().FindBin(17.25), mass_histos[0].GetXaxis().FindBin(17.75), mass_histos[0].GetXaxis().FindBin(18.25), mass_histos[0].GetXaxis().FindBin(19.25)]
    intLabel = ROOT.TLatex()
    intLabel.SetNDC()
    intLabel.SetTextSize(0.02)
    intLabel.SetTextFont(62)
    massIntegrals = [[w.Integral(lowBin[0], highBin[0]) for w in mass_histos], [w.Integral(lowBin[1], highBin[1]) for w in mass_histos], [w.Integral(lowBin[2], highBin[2]) for w in mass_histos], [w.Integral(lowBin[3], highBin[3]) for w in mass_histos]]
    
    mass_histos_scaled = [w.Clone() for w in mass_histos]
    for j in range(len(mass_histos_scaled)):
        mass_histos_scaled[j].Scale(1.0/lC)
        mass_histos_scaled[j].SetTitle("LiveCharge Normalized (1/mC) Invariant Mass of Potential Particle 0.5 MeV per bin Run %d - Cut: %s;Mass [MeV/c^{2}]; Count [1/mC]")
    massIntegrals_scaled = [[w.Integral(lowBin[0], highBin[0]) for w in mass_histos_scaled], [w.Integral(lowBin[1], highBin[1]) for w in mass_histos_scaled], [w.Integral(lowBin[2], highBin[2]) for w in mass_histos_scaled], [w.Integral(lowBin[3], highBin[3]) for w in mass_histos_scaled]]

    c = ROOT.TCanvas("c", "c", 1000, 1000)
    c.SetLeftMargin(0.15)
    latex = ROOT.TLatex()
    latex.SetNDC()                  # Use normalized coordinates (0 to 1)
    latex.SetTextAlign(22)          # 22 = centered horizontally and vertically
    latex.SetTextSize(0.04)
    latex.SetTextFont(62)
    latex.DrawLatex(0.5, 0.6, "3 Particle [e+e-e-] Candidate Document")
    latex.DrawLatex(0.5, 0.4, "Cut 0 - No Cuts")
    c.Update()
    c.Print(X_pdfName+"(")
    c.Clear()

    # Top left
    pad1 = ROOT.TPad("pad1", "", 0.00, 0.50, 0.50, 1.00)
    # Top right
    pad2 = ROOT.TPad("pad2", "", 0.50, 0.50, 1.00, 1.00)
    # Bottom (full width)
    pad3 = ROOT.TPad("pad3", "", 0.00, 0.00, 1.00, 0.50)
    for pad in (pad1, pad2, pad3):
        pad.Draw()

    # Draw in each pad
    pad1.cd()
    gPad.SetLogz(1)
    xy_histos[0].Draw("COLZ")

    pad2.cd()
    gPad.SetLogz(1)
    ETheta_histos[0].Draw("COLZ")

    pad3.cd()
    vZ_All_histos[0].Draw("HIST")

    c.Update()
    c.Print(X_pdfName)
    c.Clear()

    #{"None", "3 or More", "Timing", "Fiducial", "Cluster E", "E Sum", "Coplanarity", "X_E", "1 GEM Match", "2 GEM Match", "Vertex Z"}
    Cut_desc = ["Cut 1 - 3 or More Clusters", 
                "Cut 2 - Timing",
                "Cut 3 - Fiducial",
                "Cut 4 - Individual Cluster Energy",
                "Cut 5 - Energy Sum",
                "Cut 6 - Coplanarity between candidate and e'",
                "Cut 7 - Candidate Energy",
                "Cut 8 - 1 GEM Match Each",
                "Cut 9 - 2 GEM Match Each",
                "Cut 10 - Vertex Z Cut for All"
                ]
    detail = ["Require that there are at least 3 clusters",  
            "Require (#Deltat = t_{Max}-t_{Min}) < 16 ns.",
            "Require all clusters outside 1/2 of first",
            "Require each cluster to be more than",
            "Require that the 3-cluster energy sum is within 5 sigma",
            "Require that the candidate and e'",
            "Require that the candidate have a majority",
            "Require that all particles have",
            "Require that all particles have",
            "Require that all particles are from within"
            ]
    if fo == 1:
        detail[2] = "Require all clusters are in the first"

    detail2 = ["and form all combinations of 2 as candidates.",  
            "",
            "open layer, and not in outermost layer.",
            "70MeV and less than 0.8*E_{beam}",
            "of beam energy where #sigma = #Sigma_{i=1}^{3} #sqrt{#sigma_{E_{i}}}",
            "be coplanar within 15#circ",
            "of the beam energy (E_{X} > 0.5*E_{beam})",
            "a match on at least one of the GEMs",
            "a match on both GEM layers",
            "2m of the target based on DCA Z reconstruction."
            ]
    
    for i in range(1, len(xy_histos)):
        latex.DrawLatex(0.5, 0.6, Cut_desc[i-1])
        latex.DrawLatex(0.5, 0.4, detail[i-1])
        latex.DrawLatex(0.5, 0.3, detail2[i-1])
        c.Update()
        c.Print(X_pdfName)
        c.Clear()

        #HyCal Position Page
        c.Divide(2,2)
        latex.DrawLatex(0.5, 0.5, "HyCal Position")
        c.cd(1)
        gPad.SetLogz(1)
        xy_histos[i].Draw("COLZ")
        c.cd(2)
        gPad.SetLogz(1)
        ETheta_histos[i].Draw("COLZ")
        c.cd(3)
        xy_ProjX[i].Draw("HIST")
        c.cd(4)
        xy_ProjY[i].Draw("HIST")
        c.Print(X_pdfName)
        c.Clear()

        #Momentum Page
        c.Divide(2,2)
        latex.DrawLatex(0.5, 0.5, "Momentum - P")
        c.cd(1)
        gPad.SetLogz(1)
        sumPxVPy_histos[i].Draw("COLZ")
        c.cd(2)
        gPad.SetLogz(1)
        sumPt_histos[i].Draw("HIST")
        c.cd(3)
        PxvPy_ProjX[i].SetTitle(f"P_{{x}} Cut - {i}")
        PxvPy_ProjX[i].Draw("HIST")
        c.cd(4)
        PxvPy_ProjY[i].SetTitle(f"P_{{y}} Cut - {i}")
        PxvPy_ProjY[i].Draw("HIST")
        c.Print(X_pdfName)
        c.Clear()

        #Energy Page
        c.Divide(2,2)
        latex.DrawLatex(0.5, 0.5, "Energy")
        c.cd(1)
        sumE_histos[i].Draw("HIST")
        c.cd(2)
        minE_histos[i].Draw("HIST")
        c.cd(3)
        medE_histos[i].Draw("HIST")
        c.cd(4)
        maxE_histos[i].Draw("HIST")
        c.Print(X_pdfName)
        c.Clear()

        #Timing and Delta Phi
        c.Divide(1,2)
        latex.DrawLatex(0.5, 0.52, "Timing and #Delta#phi")
        c.cd(1)
        time_histos[i].Draw("HIST")
        c.cd(2)
        diffPhi_histos[i].Draw("HIST")
        c.Print(X_pdfName)
        c.Clear()

        #Vertex Z
        c.Divide(1,2)
        latex.DrawLatex(0.5, 0.52, "Vertex Z (DCA)")
        c.cd(1)
        vZ_All_histos[i].Draw("HIST")
        c.cd(2)
        vZ_MinE_histos[i].Draw("HIST")
        c.Print(X_pdfName)
        c.Clear()

        c.Divide(1,2)
        latex.DrawLatex(0.5, 0.52, "Vertex Z (DCA)")
        c.cd(1)
        vZ_MedE_histos[i].Draw("HIST")
        c.cd(2)
        vZ_MaxE_histos[i].Draw("HIST")
        c.Print(X_pdfName)
        c.Clear()

        #Invariant Mass - Counts
        mass_histos[i].Draw("HIST")
        intLabel.DrawLatex(0.4, 0.8, f"17 MeV #pm 0.5 = {massIntegrals[0][i]:,.1f}")
        intLabel.DrawLatex(0.4, 0.75, f"17 MeV #pm 1.0 = {massIntegrals[1][i]:,.1f}")
        intLabel.DrawLatex(0.4, 0.7, f"17 MeV #pm 1.5 = {massIntegrals[2][i]:,.1f}")
        intLabel.DrawLatex(0.4, 0.65, f"17 MeV #pm 2.5 = {massIntegrals[3][i]:,.1f}")
        c.Print(X_pdfName)
        c.Clear()

        #Invariant Mass - LiveCharge Normalized
        if(lC < 0.9):   
            mass_histos_scaled[i].Draw("HIST")
            intLabel.DrawLatex(0.4, 0.8, f"17 MeV #pm 0.5 = {massIntegrals_scaled[0][i]:,.1f}")
            intLabel.DrawLatex(0.4, 0.75, f"17 MeV #pm 1.0 = {massIntegrals_scaled[1][i]:,.1f}")
            intLabel.DrawLatex(0.4, 0.7, f"17 MeV #pm 1.5 = {massIntegrals_scaled[2][i]:,.1f}")
            intLabel.DrawLatex(0.4, 0.65, f"17 MeV #pm 2.5 = {massIntegrals_scaled[3][i]:,.1f}")
        c.Print(X_pdfName)
        c.Clear()
        
    c.Print(X_pdfName + ")")
    return 0

def printXPdf_reduced(X_histos, X_pdfName, lC, fo):
    #Separate out the histograms for the X candidates
    xy_histos = [w for w in X_histos if "HC_XY_" in w.GetName()]
    xy_ProjX = [w.ProjectionX() for w in xy_histos]
    xy_ProjY = [w.ProjectionY() for w in xy_histos]

    ETheta_histos = [w for w in X_histos if "E_theta_" in w.GetName()]
    sumPt_histos = [w for w in X_histos if "Sum_pt_" in w.GetName()]
    sumPxVPy_histos = [w for w in X_histos if "_pxVpy_" in w.GetName()]
    PxvPy_ProjX = [w.ProjectionX() for w in sumPxVPy_histos]
    PxvPy_ProjY = [w.ProjectionY() for w in sumPxVPy_histos]

    diffPhi_histos = [w for w in X_histos if "diffPhi_" in w.GetName()]
    time_histos = [w for w in X_histos if "_timing_" in w.GetName()]
    sumE_histos = [w for w in X_histos if "_sumE_" in w.GetName()]
    minE_histos = [w for w in X_histos if "_minE_" in w.GetName()]
    medE_histos = [w for w in X_histos if "_medE_" in w.GetName()]
    maxE_histos = [w for w in X_histos if "_maxE_" in w.GetName()]
    
    #Separate out the vertex Z of histograms for the X candidates
    vZ_histos = [w for w in X_histos if "_vZ_" in w.GetName()]
    vZ_All_histos = [w for w in vZ_histos if "_All_" in w.GetName()]
    vZ_MinE_histos = [w for w in vZ_histos if "_Min_E_" in w.GetName()]
    vZ_MedE_histos = [w for w in vZ_histos if "_Med_E_" in w.GetName()]
    vZ_MaxE_histos = [w for w in vZ_histos if "_Max_E_" in w.GetName()]

    mass_histos = [w for w in X_histos if "_invM_" in w.GetName()]
    lowBin = [mass_histos[0].GetXaxis().FindBin(16.75), mass_histos[0].GetXaxis().FindBin(16.25), mass_histos[0].GetXaxis().FindBin(15.75), mass_histos[0].GetXaxis().FindBin(14.75)]
    highBin = [mass_histos[0].GetXaxis().FindBin(17.25), mass_histos[0].GetXaxis().FindBin(17.75), mass_histos[0].GetXaxis().FindBin(18.25), mass_histos[0].GetXaxis().FindBin(19.25)]
    intLabel = ROOT.TLatex()
    intLabel.SetNDC()
    intLabel.SetTextSize(0.03)
    intLabel.SetTextFont(62)
    massIntegrals = [[w.Integral(lowBin[0], highBin[0]) for w in mass_histos], [w.Integral(lowBin[1], highBin[1]) for w in mass_histos], [w.Integral(lowBin[2], highBin[2]) for w in mass_histos], [w.Integral(lowBin[3], highBin[3]) for w in mass_histos]]
    
    mass_histos_scaled = [w.Clone() for w in mass_histos]
    for j in range(len(mass_histos_scaled)):
        mass_histos_scaled[j].Scale(1.0/lC)
        mass_histos_scaled[j].SetTitle("LiveCharge Normalized (1/mC) Invariant Mass of Potential Particle 0.5 MeV per bin Run %d - Cut: %s;Mass [MeV/c^{2}]; Count [1/mC]")
    massIntegrals_scaled = [[w.Integral(lowBin[0], highBin[0]) for w in mass_histos_scaled], [w.Integral(lowBin[1], highBin[1]) for w in mass_histos_scaled], [w.Integral(lowBin[2], highBin[2]) for w in mass_histos_scaled], [w.Integral(lowBin[3], highBin[3]) for w in mass_histos_scaled]]

    c = ROOT.TCanvas("c", "c", 1000, 1000)
    c.SetLeftMargin(0.15)
    latex = ROOT.TLatex()
    latex.SetNDC()                  # Use normalized coordinates (0 to 1)
    latex.SetTextAlign(22)          # 22 = centered horizontally and vertically
    latex.SetTextSize(0.04)
    latex.SetTextFont(62)
    latex.DrawLatex(0.5, 0.6, "3 Particle [e+e-e-] Candidate Document")
    latex.DrawLatex(0.5, 0.4, "Cut 0 - No Cuts")
    c.Update()
    c.Print(X_pdfName+"(")
    c.Clear()

    # Top left
    pad1 = ROOT.TPad("pad1", "", 0.00, 0.50, 0.50, 1.00)
    # Top right
    pad2 = ROOT.TPad("pad2", "", 0.50, 0.50, 1.00, 1.00)
    # Bottom (full width)
    pad3 = ROOT.TPad("pad3", "", 0.00, 0.00, 1.00, 0.50)
    for pad in (pad1, pad2, pad3):
        pad.Draw()

    # Draw in each pad
    pad1.cd()
    gPad.SetLogz(1)
    xy_histos[0].Draw("COLZ")

    pad2.cd()
    gPad.SetLogz(1)
    ETheta_histos[0].Draw("COLZ")

    pad3.cd()
    vZ_All_histos[0].Draw("HIST")

    c.Update()
    c.Print(X_pdfName)
    c.Clear()

    #{"None", "3 or More", "Timing", "Fiducial", "Cluster E", "E Sum", "Coplanarity", "X_E", "1 GEM Match", "2 GEM Match", "Vertex Z"}
    Cut_desc = ["Cut 1 - 3 or More Clusters", 
                "Cut 2 - Timing",
                "Cut 3 - Fiducial",
                "Cut 4 - Individual Cluster Energy",
                "Cut 5 - Energy Sum",
                "Cut 6 - Coplanarity between candidate and e'",
                "Cut 7 - Candidate Energy",
                "Cut 8 - 1 GEM Match Each",
                "Cut 9 - 2 GEM Match Each",
                "Cut 10 - Vertex Z Cut for All"
                ]
    detail = ["Require that there are at least 3 clusters",  
            "Require (#Deltat = t_{Max}-t_{Min}) < 16 ns.",
            "Require all clusters outside 1/2 of first",
            "Require each cluster to be more than",
            "Require that the 3-cluster energy sum is within 5 sigma",
            "Require that the candidate and e'",
            "Require that the candidate have a majority",
            "Require that all particles have",
            "Require that all particles have",
            "Require that all particles are from within"
            ]
    if fo == 0:
        detail[2] = "Require all clusters are in the first"

    detail2 = ["and form all combinations of 2 as candidates.",  
            "",
            "open layer, and not in outermost layer.",
            "70MeV and less than 0.8*E_{beam}",
            "of beam energy where #sigma = #Sigma_{i=1}^{3} #sqrt{#sigma_{E_{i}}}",
            "be coplanar within 15#circ",
            "of the beam energy (E_{X} > 0.5*E_{beam})",
            "a match on at least one of the GEMs",
            "a match on both GEM layers",
            "2m of the target based on DCA Z reconstruction."
            ]
    
    for i in range(1, len(xy_histos)):
        latex.DrawLatex(0.5, 0.6, Cut_desc[i-1])
        latex.DrawLatex(0.5, 0.4, detail[i-1])
        latex.DrawLatex(0.5, 0.3, detail2[i-1])
        c.Update()
        c.Print(X_pdfName)
        c.Clear()
        
        if i == 3:
            #HyCal Position Page
            c.Divide(2,2)
            latex.DrawLatex(0.5, 0.5, "HyCal Position")
            c.cd(1)
            gPad.SetLogz(1)
            xy_histos[i].Draw("COLZ")
            c.cd(2)
            gPad.SetLogz(1)
            ETheta_histos[i].Draw("COLZ")
            c.cd(3)
            xy_ProjX[i].Draw("HIST")
            c.cd(4)
            xy_ProjY[i].Draw("HIST")
            c.Print(X_pdfName)
            c.Clear()

        #Momentum Page
        c.Divide(2,2)
        latex.DrawLatex(0.5, 0.5, "Momentum - P")
        c.cd(1)
        gPad.SetLogz(1)
        sumPxVPy_histos[i].Draw("COLZ")
        c.cd(2)
        gPad.SetLogz(1)
        sumPt_histos[i].Draw("HIST")
        c.cd(3)
        PxvPy_ProjX[i].SetTitle(f"P_{{x}} Cut - {i}")
        PxvPy_ProjX[i].Draw("HIST")
        c.cd(4)
        PxvPy_ProjY[i].SetTitle(f"P_{{y}} Cut - {i}")
        PxvPy_ProjY[i].Draw("HIST")
        c.Print(X_pdfName)
        c.Clear()

        #Energy Page
        c.Divide(2,2)
        latex.DrawLatex(0.5, 0.5, "Energy")
        c.cd(1)
        sumE_histos[i].Draw("HIST")
        c.cd(2)
        minE_histos[i].Draw("HIST")
        c.cd(3)
        medE_histos[i].Draw("HIST")
        c.cd(4)
        maxE_histos[i].Draw("HIST")
        c.Print(X_pdfName)
        c.Clear()

        #Timing Delta Phi and 
        # Top left
        pad1 = ROOT.TPad(f"pad1_{i}", "", 0.00, 0.50, 0.50, 1.00)
        # Top right
        pad2 = ROOT.TPad(f"pad2_{i}", "", 0.50, 0.50, 1.00, 1.00)
        # Bottom (full width)
        pad3 = ROOT.TPad(f"pad3_{i}", "", 0.00, 0.00, 1.00, 0.50)
        for pad in (pad1, pad2, pad3):
            pad.Draw()
        latex.DrawLatex(0.15, 0.5, "Timing, #Delta#phi, V_{z}")
        pad1.cd()
        time_histos[i].Draw("HIST")
        pad2.cd()
        diffPhi_histos[i].Draw("HIST")
        pad3.cd()
        vZ_All_histos[i].Draw("HIST")
        c.Update()
        c.Print(X_pdfName)
        c.Clear()

        if(i >= 8):
            c.Divide(1,3)
            c.cd(1)
            vZ_MinE_histos[i].Draw("HIST")
            c.cd(2)
            vZ_MedE_histos[i].Draw("HIST")
            c.cd(3)
            vZ_MaxE_histos[i].Draw("HIST")
            c.Print(X_pdfName)
            c.Clear()


        #Invariant Mass - Counts
        c.Divide(1,2)
        c.cd(1)
        mass_histos[i].Draw("HIST")
        intLabel.DrawLatex(0.4, 0.8, f"17 MeV #pm 0.5 = {massIntegrals[0][i]:,.1f}")
        intLabel.DrawLatex(0.4, 0.75, f"17 MeV #pm 1.0 = {massIntegrals[1][i]:,.1f}")
        intLabel.DrawLatex(0.4, 0.7, f"17 MeV #pm 1.5 = {massIntegrals[2][i]:,.1f}")
        intLabel.DrawLatex(0.4, 0.65, f"17 MeV #pm 2.5 = {massIntegrals[3][i]:,.1f}")
        #Invariant Mass - LiveCharge Normalized
        if (lC < 0.9):
            c.cd(2)
            mass_histos_scaled[i].Draw("HIST")
            intLabel.DrawLatex(0.4, 0.8, f"17 MeV #pm 0.5 = {massIntegrals_scaled[0][i]:,.1f}")
            intLabel.DrawLatex(0.4, 0.75, f"17 MeV #pm 1.0 = {massIntegrals_scaled[1][i]:,.1f}")
            intLabel.DrawLatex(0.4, 0.7, f"17 MeV #pm 1.5 = {massIntegrals_scaled[2][i]:,.1f}")
            intLabel.DrawLatex(0.4, 0.65, f"17 MeV #pm 2.5 = {massIntegrals_scaled[3][i]:,.1f}")
        c.Print(X_pdfName)
        c.Clear()
        
    c.Print(X_pdfName + ")")
    return 0

def printMPdf(M_histos, M_pdfName, lC, fo):

    #Separate out the histograms for the X candidates
    xy_histos = [w for w in M_histos if "HC_XY_" in w.GetName()]
    xy_ProjX = [w.ProjectionX() for w in xy_histos]
    xy_ProjY = [w.ProjectionY() for w in xy_histos]

    ETheta_histos = [w for w in M_histos if "E_theta_" in w.GetName()]
    sumPt_histos = [w for w in M_histos if "Sum_pt_" in w.GetName()]
    sumPxVPy_histos = [w for w in M_histos if "_pxVpy_" in w.GetName()]
    PxvPy_ProjX = [w.ProjectionX() for w in sumPxVPy_histos]
    PxvPy_ProjY = [w.ProjectionY() for w in sumPxVPy_histos]

    diffPhi_histos = [w for w in M_histos if "diffPhi_" in w.GetName()]
    time_histos = [w for w in M_histos if "_timing_" in w.GetName()]
    sumE_histos = [w for w in M_histos if "_sumE_" in w.GetName()]
    
    #Separate out the vertex Z of histograms for the X candidates
    vZ_histos = [w for w in M_histos if "_vZ_" in w.GetName()]

    mass_histos = [w for w in M_histos if "_invM_" in w.GetName()]

    c = ROOT.TCanvas("c", "c", 1000, 1000)
    c.SetLeftMargin(0.15)
    latex = ROOT.TLatex()
    latex.SetNDC()                  # Use normalized coordinates (0 to 1)
    latex.SetTextAlign(22)          # 22 = centered horizontally and vertically
    latex.SetTextSize(0.04)
    latex.SetTextFont(62)
    latex.DrawLatex(0.5, 0.6, "Moller Candidate Document")
    c.Update()
    c.Print(M_pdfName+"(")
    c.Clear()


    Cut_desc = ["Cut 1 - 2 or More Clusters", 
                "Cut 2 - Timing",
                "Cut 3 - Fiducial",
                "Cut 4 - Coplanarity",
                "Cut 5 - Elasticity",
                "Cut 6 - 1 GEM Match Each",
                "Cut 7 - 2 GEM Match Each",
                "Cut 8 - Vertex Z Cut for All"
                ]
    detail = ["Require that there are at least 2 clusters.",  
            "Require (#Deltat = t_{Max}-t_{Min}) < 16 ns.",
            "Require all clusters outside 1/2 of first",
            "Require that each cluster to",
            "Require a check for elasticity",
            "Require that all particles have",
            "Require that all particles have",
            "Require that all particles are from within"
            ]
    if fo == 3 or fo == 4:
        Cut_desc[0] = "Cut 1 - Exactly 2 Clusters"
        detail[0] = "Require that there are EXACTLY 2 clusters"
    if fo == 1:
        detail[2] = "Require all clusters are in the first"
    
    detail2 = ["",  
            "",
            "open layer, and not in outermost layer.",
            "be coplanar within 10#circ",
            "",
            "a match on at least one of the GEMs",
            "a match on both GEM layers",
            "2m of the target based on DCA Z reconstruction."
            ]
    
    for i in range(0, len(xy_histos)):
        latex.DrawLatex(0.5, 0.6, Cut_desc[i])
        latex.DrawLatex(0.5, 0.4, detail[i])
        latex.DrawLatex(0.5, 0.3, detail2[i])
        c.Update()
        c.Print(M_pdfName)
        c.Clear()

        #HyCal Position Page
        c.Divide(2,2)
        latex.DrawLatex(0.5, 0.5, "HyCal Position")
        c.cd(1)
        gPad.SetLogz(1)
        xy_histos[i].Draw("COLZ")
        c.cd(2)
        gPad.SetLogz(1)
        ETheta_histos[i].Draw("COLZ")
        c.cd(3)
        xy_ProjX[i].Draw("HIST")
        c.cd(4)
        xy_ProjY[i].Draw("HIST")
        c.Print(M_pdfName)
        c.Clear()

        #Momentum Page
        c.Divide(2,2)
        latex.DrawLatex(0.5, 0.5, "Momentum - P")
        c.cd(1)
        gPad.SetLogz(1)
        sumPxVPy_histos[i].Draw("COLZ")
        c.cd(2)
        gPad.SetLogz(1)
        sumPt_histos[i].Draw("HIST")
        c.cd(3)
        PxvPy_ProjX[i].SetTitle(f"P_{{x}} Cut - {i}")
        PxvPy_ProjX[i].Draw("HIST")
        c.cd(4)
        PxvPy_ProjY[i].SetTitle(f"P_{{y}} Cut - {i}")
        PxvPy_ProjY[i].Draw("HIST")
        c.Print(M_pdfName)
        c.Clear()

        #Energy, Timing, and Delta Phi Page
        c.Divide(2,2)
        latex.DrawLatex(0.5, 0.5, "Energy, Timing, and #Delta#phi")
        c.cd(1)
        sumE_histos[i].Draw("HIST")
        c.cd(2)
        time_histos[i].Draw("HIST")
        c.cd(3)
        diffPhi_histos[i].Draw("HIST")
        c.Print(M_pdfName)
        c.Clear()

        #Vertex Z and Invariant Mass
        c.Divide(1,2)
        latex.DrawLatex(0.5, 0.52, "Vertex Z (DCA) and Invariant Mass")
        c.cd(1)
        vZ_histos[i].Draw("HIST")
        c.cd(2)
        #mass_histos[i].Fit("gaus", "Q", "", 42.0, 50.0)
        #fit = mass_histos[i].GetFunction("gaus")
        if i > 5:
            fit = ROOT.TF1("cb", "crystalball", 35.0, 50.0)
            # Parameters: 0=Constant, 1=Mean, 2=Sigma, 3=Alpha, 4=N
            fit.SetParameters(mass_histos[i].GetMaximum(), mass_histos[i].GetMean(), mass_histos[i].GetRMS(), 1.0, 2.0)
            mass_histos[i].Fit("cb", "SQ")
            gStyle.SetOptStat(11)
            gStyle.SetOptFit(10)
        mass_histos[i].Draw("HIST")
        if i > 5:
            fit.Draw("SAME")
            st = mass_histos[i].GetListOfFunctions().FindObject("stats")
            if st:
                st.SetX1NDC(0.1)
                st.SetX2NDC(0.30)
                st.SetY1NDC(0.70)
                st.SetY2NDC(0.90)
                st.Draw("SAME")
        c.Print(M_pdfName)
        c.Clear()
        
    c.Print(M_pdfName + ")")

    return 0

def printXGGPdf_reduced(X_histos, X_pdfName, lC, fo):
    #Separate out the histograms for the X candidates
    xy_histos = [w for w in X_histos if "HC_XY_" in w.GetName()]
    xy_ProjX = [w.ProjectionX() for w in xy_histos]
    xy_ProjY = [w.ProjectionY() for w in xy_histos]

    ETheta_histos = [w for w in X_histos if "E_theta_" in w.GetName()]
    sumPt_histos = [w for w in X_histos if "Sum_pt_" in w.GetName()]
    sumPxVPy_histos = [w for w in X_histos if "_pxVpy_" in w.GetName()]
    PxvPy_ProjX = [w.ProjectionX() for w in sumPxVPy_histos]
    PxvPy_ProjY = [w.ProjectionY() for w in sumPxVPy_histos]

    diffPhi_histos = [w for w in X_histos if "diffPhi_" in w.GetName()]
    time_histos = [w for w in X_histos if "_timing_" in w.GetName()]
    sumE_histos = [w for w in X_histos if "_sumE_" in w.GetName()]
    minE_histos = [w for w in X_histos if "_minE_" in w.GetName()]
    medE_histos = [w for w in X_histos if "_medE_" in w.GetName()]
    maxE_histos = [w for w in X_histos if "_maxE_" in w.GetName()]
    
    #Separate out the vertex Z of histograms for the X candidates
    vZ_histos = [w for w in X_histos if "_vZ_" in w.GetName()]
    vZ_All_histos = [w for w in vZ_histos if "_All_" in w.GetName()]
    vZ_MinE_histos = [w for w in vZ_histos if "_Min_E_" in w.GetName()]
    vZ_MedE_histos = [w for w in vZ_histos if "_Med_E_" in w.GetName()]
    vZ_MaxE_histos = [w for w in vZ_histos if "_Max_E_" in w.GetName()]

    mass_histos = [w for w in X_histos if "_invM_" in w.GetName()]
    lowBin = [mass_histos[0].GetXaxis().FindBin(16.75), mass_histos[0].GetXaxis().FindBin(16.25), mass_histos[0].GetXaxis().FindBin(15.75), mass_histos[0].GetXaxis().FindBin(14.75)]
    highBin = [mass_histos[0].GetXaxis().FindBin(17.25), mass_histos[0].GetXaxis().FindBin(17.75), mass_histos[0].GetXaxis().FindBin(18.25), mass_histos[0].GetXaxis().FindBin(19.25)]
    intLabel = ROOT.TLatex()
    intLabel.SetNDC()
    intLabel.SetTextSize(0.03)
    intLabel.SetTextFont(62)
    massIntegrals = [[w.Integral(lowBin[0], highBin[0]) for w in mass_histos], [w.Integral(lowBin[1], highBin[1]) for w in mass_histos], [w.Integral(lowBin[2], highBin[2]) for w in mass_histos], [w.Integral(lowBin[3], highBin[3]) for w in mass_histos]]
    
    mass_histos_scaled = [w.Clone() for w in mass_histos]
    for j in range(len(mass_histos_scaled)):
        mass_histos_scaled[j].Scale(1.0/lC)
        mass_histos_scaled[j].SetTitle("LiveCharge Normalized (1/mC) Invariant Mass of Potential Particle 0.5 MeV per bin Run %d - Cut: %s;Mass [MeV/c^{2}]; Count [1/mC]")
    massIntegrals_scaled = [[w.Integral(lowBin[0], highBin[0]) for w in mass_histos_scaled], [w.Integral(lowBin[1], highBin[1]) for w in mass_histos_scaled], [w.Integral(lowBin[2], highBin[2]) for w in mass_histos_scaled], [w.Integral(lowBin[3], highBin[3]) for w in mass_histos_scaled]]

    c = ROOT.TCanvas("c", "c", 1000, 1000)
    c.SetLeftMargin(0.15)
    latex = ROOT.TLatex()
    latex.SetNDC()                  # Use normalized coordinates (0 to 1)
    latex.SetTextAlign(22)          # 22 = centered horizontally and vertically
    latex.SetTextSize(0.04)
    latex.SetTextFont(62)
    latex.DrawLatex(0.5, 0.6, "3 Particle [#gamma#gammae-] Candidate Document")
    latex.DrawLatex(0.5, 0.4, "Cut 0 - No Cuts")
    c.Update()
    c.Print(X_pdfName+"(")
    c.Clear()

    # Top left
    pad1 = ROOT.TPad("pad1", "", 0.00, 0.50, 0.50, 1.00)
    # Top right
    pad2 = ROOT.TPad("pad2", "", 0.50, 0.50, 1.00, 1.00)
    # Bottom (full width)
    pad3 = ROOT.TPad("pad3", "", 0.00, 0.00, 1.00, 0.50)
    for pad in (pad1, pad2, pad3):
        pad.Draw()

    # Draw in each pad
    pad1.cd()
    gPad.SetLogz(1)
    xy_histos[0].Draw("COLZ")

    pad2.cd()
    gPad.SetLogz(1)
    ETheta_histos[0].Draw("COLZ")

    pad3.cd()
    vZ_All_histos[0].Draw("HIST")

    c.Update()
    c.Print(X_pdfName)
    c.Clear()

    #{"None", "3 or More", "Timing", "Fiducial", "Cluster E", "E Sum", "Coplanarity", "X_E", "1 GEM Match", "2 GEM Match", "Vertex Z"}
    Cut_desc = ["Cut 1 - 3 or More Clusters", 
                "Cut 2 - Timing",
                "Cut 3 - Fiducial",
                "Cut 4 - Individual Cluster Energy",
                "Cut 5 - Energy Sum",
                "Cut 6 - Coplanarity between candidate and e'",
                "Cut 7 - Candidate Energy",
                "Cut 8 - No GEM matches for decay #gamma#gamma",
                "Cut 9 - 1 GEM Match for e'",
                "Cut 10 - 2 GEM Match for e'",
                "Cut 11 - Vertex Z Cut for All"
                ]
    detail = ["Require that there are at least 3 clusters",  
            "Require (#Deltat = t_{Max}-t_{Min}) < 16 ns.",
            "Require all clusters outside 1/2 of first",
            "Require each cluster to be more than",
            "Require that the 3-cluster energy sum is within 5 sigma",
            "Require that the candidate and e'",
            "Require that the candidate have a majority",
            "Require that the decay #gamma#gamma",
            "Require that all particles have",
            "Require that all particles have",
            "Require that all particles are from within"
            ]
    if fo == 0:
        detail[2] = "Require all clusters are in the first"

    detail2 = ["and form all combinations of 2 as candidates.",  
            "",
            "open layer, and not in outermost layer.",
            "70MeV and less than 0.8*E_{beam}",
            "of beam energy where #sigma = #Sigma_{i=1}^{3} #sqrt{#sigma_{E_{i}}}",
            "be coplanar within 15#circ",
            "of the beam energy (E_{X} > 0.5*E_{beam})",
            "have no matches on any GEM",
            "a match on at least one of the GEMs",
            "a match on both GEM layers",
            "2m of the target based on DCA Z reconstruction."
            ]
    
    for i in range(1, len(xy_histos)):
        latex.DrawLatex(0.5, 0.6, Cut_desc[i-1])
        latex.DrawLatex(0.5, 0.4, detail[i-1])
        latex.DrawLatex(0.5, 0.3, detail2[i-1])
        c.Update()
        c.Print(X_pdfName)
        c.Clear()
        
        if True or i == 3:
            #HyCal Position Page
            c.Divide(2,2)
            latex.DrawLatex(0.5, 0.5, "HyCal Position")
            c.cd(1)
            gPad.SetLogz(1)
            xy_histos[i].Draw("COLZ")
            c.cd(2)
            gPad.SetLogz(1)
            ETheta_histos[i].Draw("COLZ")
            c.cd(3)
            xy_ProjX[i].Draw("HIST")
            c.cd(4)
            xy_ProjY[i].Draw("HIST")
            c.Print(X_pdfName)
            c.Clear()

        #Momentum Page
        c.Divide(2,2)
        latex.DrawLatex(0.5, 0.5, "Momentum - P")
        c.cd(1)
        gPad.SetLogz(1)
        sumPxVPy_histos[i].Draw("COLZ")
        c.cd(2)
        gPad.SetLogz(1)
        sumPt_histos[i].Draw("HIST")
        c.cd(3)
        PxvPy_ProjX[i].SetTitle(f"P_{{x}} Cut - {i}")
        PxvPy_ProjX[i].Draw("HIST")
        c.cd(4)
        PxvPy_ProjY[i].SetTitle(f"P_{{y}} Cut - {i}")
        PxvPy_ProjY[i].Draw("HIST")
        c.Print(X_pdfName)
        c.Clear()

        #Energy Page
        c.Divide(2,2)
        latex.DrawLatex(0.5, 0.5, "Energy")
        c.cd(1)
        sumE_histos[i].Draw("HIST")
        c.cd(2)
        minE_histos[i].Draw("HIST")
        c.cd(3)
        medE_histos[i].Draw("HIST")
        c.cd(4)
        maxE_histos[i].Draw("HIST")
        c.Print(X_pdfName)
        c.Clear()

        #Timing Delta Phi and 
        # Top left
        pad1 = ROOT.TPad(f"pad1_{i}", "", 0.00, 0.50, 0.50, 1.00)
        # Top right
        pad2 = ROOT.TPad(f"pad2_{i}", "", 0.50, 0.50, 1.00, 1.00)
        # Bottom (full width)
        pad3 = ROOT.TPad(f"pad3_{i}", "", 0.00, 0.00, 1.00, 0.50)
        for pad in (pad1, pad2, pad3):
            pad.Draw()
        latex.DrawLatex(0.15, 0.5, "Timing, #Delta#phi, V_{z}")
        pad1.cd()
        time_histos[i].Draw("HIST")
        pad2.cd()
        diffPhi_histos[i].Draw("HIST")
        pad3.cd()
        vZ_All_histos[i].Draw("HIST")
        c.Update()
        c.Print(X_pdfName)
        c.Clear()

        if(i >= 7):
            c.Divide(1,3)
            c.cd(1)
            vZ_MinE_histos[i].Draw("HIST")
            c.cd(2)
            vZ_MedE_histos[i].Draw("HIST")
            c.cd(3)
            vZ_MaxE_histos[i].Draw("HIST")
            c.Print(X_pdfName)
            c.Clear()


        #Invariant Mass - Counts
        c.Divide(1,2)
        c.cd(1)
        mass_histos[i].Draw("HIST")
        intLabel.DrawLatex(0.4, 0.8, f"17 MeV #pm 0.5 = {massIntegrals[0][i]:,.1f}")
        intLabel.DrawLatex(0.4, 0.75, f"17 MeV #pm 1.0 = {massIntegrals[1][i]:,.1f}")
        intLabel.DrawLatex(0.4, 0.7, f"17 MeV #pm 1.5 = {massIntegrals[2][i]:,.1f}")
        intLabel.DrawLatex(0.4, 0.65, f"17 MeV #pm 2.5 = {massIntegrals[3][i]:,.1f}")
        #Invariant Mass - LiveCharge Normalized
        if (lC < 0.9):
            c.cd(2)
            mass_histos_scaled[i].Draw("HIST")
            intLabel.DrawLatex(0.4, 0.8, f"17 MeV #pm 0.5 = {massIntegrals_scaled[0][i]:,.1f}")
            intLabel.DrawLatex(0.4, 0.75, f"17 MeV #pm 1.0 = {massIntegrals_scaled[1][i]:,.1f}")
            intLabel.DrawLatex(0.4, 0.7, f"17 MeV #pm 1.5 = {massIntegrals_scaled[2][i]:,.1f}")
            intLabel.DrawLatex(0.4, 0.65, f"17 MeV #pm 2.5 = {massIntegrals_scaled[3][i]:,.1f}")
        c.Print(X_pdfName)
        c.Clear()
        
    c.Print(X_pdfName + ")")
    return 0

# Main execution
if __name__ == "__main__":

    #if(sys.argc == 2):
        #print("Please input only the root file to print.")
    file_path = sys.argv[1]

    lC = 1
    if len(sys.argv) >= 3:
        live_chargeFile = sys.argv[2]
        fp = Path(live_chargeFile)
        if fp.is_file():
            with open(live_chargeFile) as f:
                liveChargeData = json.load(f)
            lC = liveChargeData.get("value_nC")/1000000.0 #get the livecharge in mC
    print()
    print(lC)
    print()
    fiducial_option = 0
    if len(sys.argv) >= 4:
        fiducial_option = int(sys.argv[3])
    
    pdfName = file_path
    pdfName.replace(".root", "")
    pdfBase = pdfName[:pdfName.rfind(".")]
    X_pdfName = pdfBase + "_3Particle.pdf"
    M_pdfName = pdfBase + "_Moller.pdf"
    
    X_pdfName_red = pdfBase + "_3Particle_reduced.pdf"
    Xgg_pdfName_red = pdfBase + "_3Particle_Gamma_reduced.pdf"
    
    # Open the ROOT file in read mode
    root_file = ROOT.TFile.Open(file_path, "READ")
    if not root_file or root_file.IsZombie():
        print(f"Error: Could not open file {file_path}")
        exit(1)
        
    print(f"Scanning {file_path} for histograms...")
    all_histograms = extract_histograms(root_file)

    # Close the file safely
    root_file.Close()
    
    print(f"\nSuccessfully read {len(all_histograms)} histograms.")

    for name, h in all_histograms.items():
        h.GetYaxis().SetTitleOffset(1.5)
        #title = h.GetTitle()
        #newTitle = title[:title.rfind('5')] + " & 25667" + title[title.rfind('5')+1:]
        #h.SetTitle(newTitle)

    X_histos = [hist for name, hist in all_histograms.items() if "_X_" in name]
    printXPdf(X_histos, X_pdfName, lC, fiducial_option)
    printXPdf_reduced(X_histos, X_pdfName_red, lC, fiducial_option)

    Xgg_histos = [hist for name, hist in all_histograms.items() if "_Xgg_" in name]
    printXGGPdf_reduced(Xgg_histos, Xgg_pdfName_red, lC, fiducial_option)
    
    
    M_histos = [hist for name, hist in all_histograms.items() if "_M_" in name]
    printMPdf(M_histos, M_pdfName, lC, fiducial_option)


