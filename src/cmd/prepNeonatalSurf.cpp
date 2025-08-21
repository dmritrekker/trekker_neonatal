#include "cmd.h"

using namespace NIBR;

// COMBINED   0
// L_WM       1
// R_WM       2
// L_GM       3
// R_GM       4
// L_SUB      5
// R_SUB      6
// CSF        7
// L_VDC      8
// R_VDC      9
// CER_WM    10
// CER_GM    11
// BS        12
// SPINE     13
// CRAN      14

namespace CMDARGS_PREP_NEONATAL_SURF {
    std::string ibeatFolder = "";
    std::string mcribsFolder= "";
    std::string outputFolder= "";

    int numberOfThreads     = 0;
    std::string verbose     = "info";
    bool force              = false;
}

using namespace CMDARGS_PREP_NEONATAL_SURF;

Surface fsAseg2Surf(Image<int>& asegImg, int asegLabel, float meanFaceArea = 0)
{
    Surface surf = label2surface(asegImg,asegLabel,((meanFaceArea==0) ? 0.25 : meanFaceArea));
    std::vector<int> labelField(surf.nv,asegLabel);
    surf.fields.push_back(surf.makeVertField("label",labelField));
    return surf;
}

void run_prep_neonatal_surf()
{

    parseCommon(numberOfThreads,verbose);

    Image<int> asegImg(mcribsFolder + "/mri/aparc+aseg.mgz");
    asegImg.read();

    // Option 1: First combine labels and generate surface
    Image<int> cerImg;
    cerImg.createFromTemplate(asegImg,true);
    for (int n = 0; n < asegImg.numel; n++) {
        cerImg.data[n] = (asegImg.data[n]==91 || asegImg.data[n]==93) ? 1 : 0; 
    }
    Surface cer = fsAseg2Surf(cerImg,1);

    Image<int> bstImg;
    bstImg.createFromTemplate(asegImg, true);
    for (int n = 0; n < asegImg.numel; n++) {
        bstImg.data[n] = (asegImg.data[n] == 170) ? 1 : 0;
    }
    Surface bst = fsAseg2Surf(bstImg, 1);

    Image<int> sgmImg;
    sgmImg.createFromTemplate(asegImg, true);
    for (int n = 0; n < asegImg.numel; n++) {
        sgmImg.data[n] = (
            asegImg.data[n]==9 || asegImg.data[n]==11 || 
            asegImg.data[n]==13 || asegImg.data[n]==28 ||
            asegImg.data[n]==48 || asegImg.data[n]==50 ||
            asegImg.data[n]==52 || asegImg.data[n]==60 ) ? 1 : 0; 
    }
    Surface sgm = fsAseg2Surf(sgmImg, 1);

    Image<int> csfImg;
    csfImg.createFromTemplate(asegImg, true);
    for (int n = 0; n < asegImg.numel; n++) {
        csfImg.data[n] = (
            asegImg.data[n]==4 || asegImg.data[n]==14 || 
            asegImg.data[n]==15 || asegImg.data[n]==43 ) ? 1 : 0; 
    }
    Surface csf = fsAseg2Surf(csfImg, 1);

    /*  Option 2: First generate surfaces then combine
    Surface l_cer = fsAseg2Surf(cerImg,91);
    Surface r_cer = fsAseg2Surf(cerImg,93);
    Surface cer   = surfMerge(r_cer,l_cer);
    */

    makeFolder(outputFolder);

    cer.write(outputFolder + "/cer.vtk" );
    bst.write(outputFolder + "/bst.vtk" );
    sgm.write(outputFolder + "/sgm.vtk" );
    csf.write(outputFolder + "/csf.vtk" );

    // Read the original ribbon data
    Image<int> ribbonOrig(mcribsFolder + "/mri/ribbon.mgz");
    ribbonOrig.read();

    // Threshold for faster interplolation
    Image<int8_t> ribbon;
    ribbon.createFromTemplate(ribbonOrig,true);
    for (int n = 0; n < ribbonOrig.numel; n++) {
        ribbon.data[n] = (ribbonOrig.data[n]==3 || ribbonOrig.data[n]==42) ? 1 : 0; 
    }
    imgDilate(ribbon);
    imgDilate(ribbon);
    ribbon.setInterpolationMethod(NEAREST);
    // ribbon.write(outputFolder + "/ribbonThresh.nii.gz");

    /*
    Surface r_gm_closed = Surface(ibeatFolder + "/T2-iso-skullstripped-rmcere-tissue.rh.OuterSurf.PhysicalSpace.vtk");
    r_gm_closed.readMesh();
    Surface l_gm_closed = Surface(ibeatFolder + "/T2-iso-skullstripped-rmcere-tissue.lh.OuterSurf.PhysicalSpace.vtk");
    l_gm_closed.readMesh();

    Surface cortex_gm = surfMerge(r_gm_closed,l_gm_closed);
    cortex_gm.write(outputFolder + "/cortex_gm.vtk");

    Image<float> edt_gm;
    edt_gm.createFromTemplate(ribbonOrig,true);
    mapSurface2Image(&cortex_gm, &edt_gm, 0, NULL, NULL, EDT);
    edt_gm.write(outputFolder + "/edt_gm.nii.gz");



    Surface r_wm_closed = Surface(ibeatFolder + "/T2-iso-skullstripped-rmcere-tissue.rh.InnerSurf.PhysicalSpace.vtk");
    r_wm_closed.readMesh();
    Surface r_wm_closed_inf = surfMoveVerticesAlongNormal(r_wm_closed, 2);

    Surface l_wm_closed = Surface(ibeatFolder + "/T2-iso-skullstripped-rmcere-tissue.lh.InnerSurf.PhysicalSpace.vtk");
    l_wm_closed.readMesh();
    Surface l_wm_closed_inf = surfMoveVerticesAlongNormal(l_wm_closed, 2);

    Surface cortex_wm = surfMerge(r_wm_closed_inf,l_wm_closed_inf);
    cortex_wm.write(outputFolder + "/cortex_wm.vtk");

    Image<float> edt_wm;
    edt_wm.createFromTemplate(ribbonOrig,true);
    mapSurface2Image(&cortex_wm, &edt_wm, 0, NULL, NULL, EDT);
    edt_wm.write(outputFolder + "/edt_wm.nii.gz");


    for (int n = 0; n < ribbonOrig.numel; n++) {
        ribbon.data[n] = (edt_wm.data[n] < 0 && edt_gm.data[n] > 0) ? 1 : 0; 
    }
    ribbon.write(outputFolder + "/ribbonThreshSurf.nii.gz");
    */



    
    // ------
    // Create neighboring spheres in physical space
    // We will create 3 shells, with 1 mm, 2 mm, and 3 mm radii from the center, i.e. vertex
    std::vector<std::array<float,3>> nei;
    nei.reserve(180);
    for (int n = 0; n < 20; n ++) {
        std::array<float,3> p;
        p[0] = FULLSPHERE20[n][0] * 1;
        p[1] = FULLSPHERE20[n][1] * 1;
        p[2] = FULLSPHERE20[n][2] * 1;
        nei.push_back(p);
    }
    for (int n = 0; n < 60; n ++) {
        std::array<float,3> p;
        p[0] = FULLSPHERE60[n][0] * 2;
        p[1] = FULLSPHERE60[n][1] * 2;
        p[2] = FULLSPHERE60[n][2] * 2;
        nei.push_back(p);
    }
    for (int n = 0; n < 100; n ++) {
        std::array<float,3> p;
        p[0] = FULLSPHERE100[n][0] * 4;
        p[1] = FULLSPHERE100[n][1] * 4;
        p[2] = FULLSPHERE100[n][2] * 4;
        nei.push_back(p);
    }
    // ------
    

    
    // WM surface
    Surface l_wm_closed = Surface(ibeatFolder + "/T2-iso-skullstripped-rmcere-tissue.lh.InnerSurf.PhysicalSpace.vtk");
    l_wm_closed.readMesh();
    l_wm_closed.write(outputFolder + "/l_wm.vtk");
    
    Surface r_wm_closed = Surface(ibeatFolder + "/T2-iso-skullstripped-rmcere-tissue.rh.InnerSurf.PhysicalSpace.vtk");
    r_wm_closed.readMesh();
    r_wm_closed.write(outputFolder + "/r_wm.vtk");


    // Create GM ribbons

    Surface l_gm_closed = Surface(ibeatFolder + "/T2-iso-skullstripped-rmcere-tissue.lh.OuterSurf.PhysicalSpace.vtk");
    l_gm_closed.readMesh();

    // l_wm_closed.calcNormalsOfVertices();
    // l_wm_closed.write(outputFolder + "/l_wm_closed.vtk");

    std::vector<int> l_midLine;
    l_midLine.reserve(l_wm_closed.nv);

    // std::vector<int> midIdx;

    for (int n = 0; n < l_wm_closed.nv; n++) {

        bool isMidLine = true;
        
        // for (float d = -3; d < 3; d += 0.1) {
        //     float r[3];
        //     r[0] = l_wm_closed.vertices[n][0] + l_wm_closed.normalsOfVertices[n][0] * d;
        //     r[1] = l_wm_closed.vertices[n][1] + l_wm_closed.normalsOfVertices[n][1] * d;
        //     r[2] = l_wm_closed.vertices[n][2] + l_wm_closed.normalsOfVertices[n][2] * d;
        //     float val = ribbon(r);
        //     if (val > EPS6) {
        //         isMidLine = false;
        //         break;
        //     }
        // }

        for (const auto& p : nei) {
            float r[3];
            r[0] = p[0] + l_wm_closed.vertices[n][0];
            r[1] = p[1] + l_wm_closed.vertices[n][1];
            r[2] = p[2] + l_wm_closed.vertices[n][2];
            float val = ribbon(r);
            if (val > EPS6) {
                isMidLine = false;
                break;
            }
        }

        if (isMidLine) {
            for (const auto& p : nei) {
                float r[3];
                r[0] = p[0] + l_gm_closed.vertices[n][0];
                r[1] = p[1] + l_gm_closed.vertices[n][1];
                r[2] = p[2] + l_gm_closed.vertices[n][2];
                float val = ribbon(r);
                if (val > EPS6) {
                    isMidLine = false;
                    break;
                }
            }
        }

        // if (isMidLine) {
        //     midIdx.push_back(n);
        // }
        l_midLine.push_back(int(isMidLine));
    }

    SurfaceField l_midLineMask = l_wm_closed.makeVertField("midLineMask", l_midLine);

    Surface l_wm_open;
    removeVertices(&l_wm_open,&l_wm_closed, &l_midLineMask);
    // l_wm_open.write(outputFolder + "/l_wm_open.vtk");

    Surface l_gm_open;
    removeVertices(&l_gm_open,&l_gm_closed, &l_midLineMask);
    // l_gm_open.write(outputFolder + "/l_gm_open.vtk");

    Surface l_gm_ribbon = surfGlueBoundaries(l_wm_open, l_gm_open);
    l_gm_ribbon.write(outputFolder + "/l_gm.vtk");


    //
    Surface r_gm_closed = Surface(ibeatFolder + "/T2-iso-skullstripped-rmcere-tissue.rh.OuterSurf.PhysicalSpace.vtk");
    r_gm_closed.readMesh();

    std::vector<int> r_midLine;
    r_midLine.reserve(r_wm_closed.nv);

    for (int n = 0; n < r_wm_closed.nv; n++) {

        bool isMidLine = true;

        for (const auto& p : nei) {
            float r[3];
            r[0] = p[0] + r_wm_closed.vertices[n][0];
            r[1] = p[1] + r_wm_closed.vertices[n][1];
            r[2] = p[2] + r_wm_closed.vertices[n][2];
            float val = ribbon(r);
            if (val > EPS6) {
                isMidLine = false;
                break;
            }
        }

        if (isMidLine) {
            for (const auto& p : nei) {
                float r[3];
                r[0] = p[0] + r_gm_closed.vertices[n][0];
                r[1] = p[1] + r_gm_closed.vertices[n][1];
                r[2] = p[2] + r_gm_closed.vertices[n][2];
                float val = ribbon(r);
                if (val > EPS6) {
                    isMidLine = false;
                    break;
                }
            }
        }

        r_midLine.push_back(int(isMidLine));
    }

    SurfaceField r_midLineMask = r_wm_closed.makeVertField("midLineMask", r_midLine);

    Surface r_wm_open;
    removeVertices(&r_wm_open,&r_wm_closed, &r_midLineMask);
    // r_wm_open.write(outputFolder + "/r_wm_open.vtk");

    Surface r_gm_open;
    removeVertices(&r_gm_open,&r_gm_closed, &r_midLineMask);
    // r_gm_open.write(outputFolder + "/r_gm_open.vtk");

    Surface r_gm_ribbon = surfGlueBoundaries(r_wm_open, r_gm_open);
    r_gm_ribbon.write(outputFolder + "/r_gm.vtk");


    // For tractography

    Surface wm = surfMerge(r_wm_closed, l_wm_closed);
    Surface gm = surfMerge(r_gm_ribbon, l_gm_ribbon);
    
    Surface seed = surfMerge(wm, bst);
    seed = surfMerge(seed, cer);
    
    Surface discard_seed = surfMerge(gm, sgm);
    Surface req_end_inside = surfMerge(discard_seed, cer);
    req_end_inside = surfMerge(req_end_inside, bst);

    wm.write(outputFolder + "/wm.vtk");
    gm.write(outputFolder + "/gm.vtk");

    seed.write(outputFolder + "/seed.vtk");
    discard_seed.write(outputFolder + "/discard_seed.vtk");
    req_end_inside.write(outputFolder + "/req_end_inside.vtk");

    return;

}

void prepNeonatalSurf(CLI::App* app)
{
    app->description("creates surfaces for neonatal whole-brain tractography");
    
    app->add_option("<iBEAT_folder>",        ibeatFolder, "Path to iBEAT folder")
        ->required()
        ->check(CLI::ExistingDirectory);
    
    app->add_option("<MCRIBS_folder>",       mcribsFolder, "Path to MCRIBS folder")
        ->required()
        ->check(CLI::ExistingDirectory);

    app->add_option("<output_folder>",       outputFolder, "Path to output folder")
        ->required();

    app->add_option("--numberOfThreads, -n", numberOfThreads,    "Number of threads.")->check(CLI::Range(0, INT32_MAX));
    app->add_option("--verbose, -v",         verbose,            "Verbose level. Options are \"quite\",\"fatal\",\"error\",\"warn\",\"info\" and \"debug\". Default=info");
    app->add_flag("--force, -f",             force,              "Force overwriting of existing file");
    
    app->callback(run_prep_neonatal_surf);

}
