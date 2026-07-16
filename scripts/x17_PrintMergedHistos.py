import ROOT
import sys

from ROOT import kBlack
from ROOT import kBlue
from ROOT import kMagenta
from ROOT import kPink
from ROOT import kRed
from ROOT import kGreen
from ROOT import kCyan
from ROOT import kSpring
from ROOT import gPad

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

# Main execution
if __name__ == "__main__":

    #if(sys.argc == 2):
        #print("Please input only the root file to print.")
    file_path = sys.argv[1]

    pdfName = file_path
    pdfName.replace(".root", "")
    X_pdfName = "../outfiles/"+pdfName[:pdfName.rfind("/")]+"_X.pdf"
    M_pdfName = "../outfiles/"+pdfName[:pdfName.rfind("/")]+"_Moller.pdf"
    
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

    X_histos = [hist for name, hist in all_histograms.items() if "_X_" in name]

    #Separate out the histograms for the X candidates
    xy_histos = [w for w in X_histos if "HC_XY_" in w.GetName()]
    ETheta_histos = [w for w in X_histos if "E_theta_" in w.GetName()]
    sumPt_histos = [w for w in X_histos if "Sum_pt_" in w.GetName()]
    sumPxVPy_histos = [w for w in X_histos if "_pxVpy_" in w.GetName()]
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

    c = ROOT.TCanvas("c", "c", 1000, 1000)
    latex = ROOT.TLatex()
    latex.SetNDC()                  # Use normalized coordinates (0 to 1)
    latex.SetTextAlign(22)          # 22 = centered horizontally and vertically
    latex.SetTextSize(0.05)
    latex.SetTextFont(62)
    latex.DrawLatex(0.5, 0.5, "3 Particle Candidate Document - Cut 0: No Cuts")
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
                "Cut 3 - Fiducial"
                ]
    detail = ["",  
            "Require (#Deltat = t_{Max}-t_{Min}) < 16 ns",
            "Require all clusters outside 1/2 of first open layer, and not in outermost layer."
            ]
    for i in range(1, len(xy_histos)):
        
    
    
    M_histos = [hist for name, hist in all_histograms.items() if "_M_" in name]
    
    
    
    c = ROOT.TCanvas("c1", "Merged Canvas", 1000, 1000)
    legend = ROOT.TLegend(0.1,0.8,0.4,0.9)
    bitThr = ["650", "750", "950", "1300", "1100", "1400", "1500"]
    Bit_Color = [kBlue, kRed, kGreen,kCyan, kMagenta, kSpring, 28]


