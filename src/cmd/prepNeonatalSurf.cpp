#include "cmd.h"

using namespace NIBR;

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

    // Create the output folder
    makeFolder(outputFolder);

    // Read aparc+aseg.mgz
    Image<int> asegImg(mcribsFolder + "/mri/aparc+aseg.mgz");
    asegImg.read();

    // Read volumetric labels
    Image<int> labelsImg(mcribsFolder + "/mri/mcrib_labels.nii.gz");
    labelsImg.read();

    // Read ribbon.mgz
    Image<int> ribbon(mcribsFolder + "/mri/ribbon.mgz");
    ribbon.read();

    // Extract the GM ribbon only
    Image<int8_t> gm_ribbon;
    gm_ribbon.createFromTemplate(ribbon,true);
    for (int n = 0; n < ribbon.numel; n++) {
        gm_ribbon.data[n]  = (ribbon.data[n]==3 || ribbon.data[n]==42) ? 1 : 0;
    }
    imgDilate(gm_ribbon);
    imgDilate(gm_ribbon);
    gm_ribbon.setInterpolationMethod(NEAREST);
    // gm_ribbon.write(outputFolder + "/ribbonThresh.nii.gz");

    // Generate and save a tight background mask
    disp(MSG_INFO,"Generating a tight background image");
    Image<int8_t> background;
    background.createFromTemplate(ribbon,true);
    for (int n = 0; n < background.numel; n++) {
        background.data[n]  = 
            (   (ribbon.data[n] > 0) || 
                (asegImg.data[n]==91 || asegImg.data[n]==93) ||
                (asegImg.data[n]==170) )
                ? 0 : 1;
    }
    background.write(outputFolder + "/background.nii.gz" );
    disp(MSG_INFO,"background saved.");
    

    // For CER, first combine labels and generate surface
    disp(MSG_INFO,"Generating cer surface");
    Image<int> cerImg;
    cerImg.createFromTemplate(asegImg,true);
    for (int n = 0; n < asegImg.numel; n++) {
        cerImg.data[n] = (asegImg.data[n]==91 || asegImg.data[n]==93) ? 1 : 0; 
    }
    Surface cer = fsAseg2Surf(cerImg,1);
    cer.write(outputFolder + "/cer.vtk" );
    disp(MSG_INFO,"cer saved.");

    // BST surface
    disp(MSG_INFO,"Generating bst surface");
    Surface bst   = fsAseg2Surf(asegImg, 170);
    bst.write(outputFolder + "/bst.vtk" );
    disp(MSG_INFO,"bst saved.");

    // For SGM, first generate the surfaces and then combine them
    disp(MSG_INFO,"Generating sgm surface");

    Surface sgm_left = fsAseg2Surf(labelsImg, 9);         disp(MSG_INFO,"L THALAMUS done");
    Surface tmp      = fsAseg2Surf(labelsImg, 11);        sgm_left = surfMerge(sgm_left, tmp);     disp(MSG_INFO,"L CAUDATE done");
    tmp              = fsAseg2Surf(labelsImg, 13);        sgm_left = surfMerge(sgm_left, tmp);     disp(MSG_INFO,"L PALLIDUM done");
    tmp              = fsAseg2Surf(labelsImg, 12);        sgm_left = surfMerge(sgm_left, tmp);     disp(MSG_INFO,"L PUTAMEN done");
    tmp              = fsAseg2Surf(labelsImg, 17);        sgm_left = surfMerge(sgm_left, tmp);     disp(MSG_INFO,"L HIPPOCAMPUS done");
    tmp              = fsAseg2Surf(labelsImg, 18);        sgm_left = surfMerge(sgm_left, tmp);     disp(MSG_INFO,"L AMYGDALA done");
    tmp              = fsAseg2Surf(labelsImg, 26);        sgm_left = surfMerge(sgm_left, tmp);     disp(MSG_INFO,"L ACCUMBENS done");
    Surface l_vdc    = fsAseg2Surf(asegImg, 28);        sgm_left = surfMerge(sgm_left, l_vdc);   disp(MSG_INFO,"L VENTRAL DIENCEPHALON done");

    sgm_left.write(outputFolder + "/sgm_left.vtk");
    disp(MSG_INFO,"sgm_left saved.");

    Surface sgm_right = fsAseg2Surf(labelsImg, 48);       disp(MSG_INFO,"R THALAMUS done");
    tmp               = fsAseg2Surf(labelsImg, 50);       sgm_right = surfMerge(sgm_right, tmp);   disp(MSG_INFO,"R CAUDATE done");
    tmp               = fsAseg2Surf(labelsImg, 52);       sgm_right = surfMerge(sgm_right, tmp);   disp(MSG_INFO,"R PALLIDUM done");
    tmp               = fsAseg2Surf(labelsImg, 51);       sgm_right = surfMerge(sgm_right, tmp);   disp(MSG_INFO,"R PUTAMEN done");
    tmp               = fsAseg2Surf(labelsImg, 53);       sgm_right = surfMerge(sgm_right, tmp);   disp(MSG_INFO,"R HIPPOCAMPUS done");
    tmp               = fsAseg2Surf(labelsImg, 54);       sgm_right = surfMerge(sgm_right, tmp);   disp(MSG_INFO,"R AMYGDALA done");
    tmp               = fsAseg2Surf(labelsImg, 58);       sgm_right = surfMerge(sgm_right, tmp);   disp(MSG_INFO,"R ACCUMBENS done");
    Surface r_vdc     = fsAseg2Surf(asegImg, 60);       sgm_right = surfMerge(sgm_right, r_vdc); disp(MSG_INFO,"R VENTRAL DIENCEPHALON done");

    sgm_right.write(outputFolder + "/sgm_right.vtk");
    disp(MSG_INFO,"sgm_right saved.");

    Surface vdc = surfMerge(l_vdc, r_vdc);

    Surface sgm = surfMerge(sgm_left, sgm_right);
    sgm.write(outputFolder + "/sgm.vtk");
    disp(MSG_INFO,"sgm saved.");

    // For CSF, we will first combine the labels, then generate the surface, which will not be a single connected component
    disp(MSG_INFO,"Generating csf surface");
    Surface csf;

    Image<float> csfImg;
    csfImg.createFromTemplate(asegImg, true);
    for (int n = 0; n < asegImg.numel; n++) {
        csfImg.data[n] = (
              asegImg.data[n]==4  ||                            // L LATERAL VENTRICLE
              asegImg.data[n]==43 ||                            // R LATERAL VENTRICLE
              asegImg.data[n]==14 ||                            // 3rd VENTRICLE
              asegImg.data[n]==15 ||                            // 4rd VENTRICLE
              (background.data[n] == 0 && asegImg.data[n]==24)  // CSF but that is only inside the brain mask
              ) ? 1.0f : 0.0f; 
    }
    
    for (int n = 0; n < 2; n++) {
        imgDilate(csfImg,NIBR::CONN6);
        imgErode (csfImg,NIBR::CONN6);
    }

    if (!isosurface(&csfImg, 0.5, &csf)) {
        disp(MSG_ERROR, "Failed to generate surface");
        return;
    }

    if (csf.nv > 0) csf = surfSmooth(csf,2);
    float meanFaceArea = 0.25;
    if (csf.nv > 0) csf.calcArea();
    if (csf.nv > 0) csf = surfRemesh(csf,csf.area/meanFaceArea*0.5f,1,0);
    
    csf.write(outputFolder + "/csf.vtk" );
    disp(MSG_INFO,"csf saved.");
    
    // Read WM and GM surfaces
    Surface l_wm_closed = Surface(ibeatFolder + "/T2-iso-skullstripped-rmcere-tissue.lh.InnerSurf.PhysicalSpace.vtk"); l_wm_closed.readMesh();
    Surface l_gm_closed = Surface(ibeatFolder + "/T2-iso-skullstripped-rmcere-tissue.lh.OuterSurf.PhysicalSpace.vtk"); l_gm_closed.readMesh();
    Surface r_wm_closed = Surface(ibeatFolder + "/T2-iso-skullstripped-rmcere-tissue.rh.InnerSurf.PhysicalSpace.vtk"); r_wm_closed.readMesh();
    Surface r_gm_closed = Surface(ibeatFolder + "/T2-iso-skullstripped-rmcere-tissue.rh.OuterSurf.PhysicalSpace.vtk"); r_gm_closed.readMesh();

    l_wm_closed.write(outputFolder + "/l_wm.vtk");
    r_wm_closed.write(outputFolder + "/r_wm.vtk");

    // auto prepSurf = [&] (Surface& surfToPrep) {
    //     surfToPrep.enablePointCheck(2.0f);
    //     surfToPrep.prepIglAABBTree();
    //     if (surfToPrep.centersOfFaces==NULL) surfToPrep.calcCentersOfFaces();
    //     if (surfToPrep.normalsOfFaces==NULL) surfToPrep.calcNormalsOfFaces();
    //     if (surfToPrep.triangleEdge1 ==NULL) surfToPrep.calcTriangleVectors();
    // };

    // prepSurf(l_wm_closed);
    // prepSurf(l_gm_closed);
    // prepSurf(r_wm_closed);
    // prepSurf(r_gm_closed);
    // prepSurf(bst);
    // prepSurf(sgm);

    Image<float> edt_bst;
    edt_bst.createFromTemplate(ribbon,true);
    mapSurface2Image(&bst, &edt_bst, 0, NULL, NULL, EDT);

    Image<float> edt_vdc;
    edt_vdc.createFromTemplate(ribbon,true);
    mapSurface2Image(&vdc, &edt_vdc, 0, NULL, NULL, EDT);


    auto genGMribbon = [&] (Surface ipsi_wm_closed, Surface ipsi_gm_closed, Surface contra_wm_closed) -> Surface {

        // Calculate Euclidean distance to the contra lateral wm
        Image<float> edt_contra_wm;
        edt_contra_wm.createFromTemplate(ribbon,true);
        mapSurface2Image(&contra_wm_closed, &edt_contra_wm, 0, NULL, NULL, EDT);

        // Compute distance for all vertices
        std::vector<int> midLine;
        midLine.reserve(ipsi_wm_closed.nv);

        for (int n = 0; n < ipsi_wm_closed.nv; n++) {

            // // float d_ipsi_wm_2_contra_wm = std::fabs(edt_contra_wm(ipsi_wm_closed.vertices[n]));
            // float d_ipsi_wm_2_contra_wm = contra_wm_closed.squaredDistToPoint(ipsi_wm_closed.vertices[n]);
            // if (d_ipsi_wm_2_contra_wm < 1.0f) {
            //     midLine.push_back(true);
            //     continue;
            // }

            // // float d_ipsi_gm_2_bst = std::fabs(edt_bst(ipsi_gm_closed.vertices[n]));
            // float d_ipsi_gm_2_bst = bst.squaredDistToPoint(ipsi_gm_closed.vertices[n]);
            // if (d_ipsi_gm_2_bst < 1.0f) {
            //     midLine.push_back(true);
            //     continue;
            // }

            // // float d_ipsi_gm_2_sgm = std::fabs(edt_sgm(ipsi_gm_closed.vertices[n]));
            // float d_ipsi_gm_2_sgm = sgm.squaredDistToPoint(ipsi_gm_closed.vertices[n]);
            // if (d_ipsi_gm_2_sgm < 1.0f) {
            //     midLine.push_back(true);
            //     continue;
            // }
            
            float d_ipsi_wm_2_contra_wm = std::fabs(edt_contra_wm(ipsi_wm_closed.vertices[n]));
            if (d_ipsi_wm_2_contra_wm < 1.0f) {
                midLine.push_back(true);
                continue;
            }

            float d_ipsi_gm_2_bst = std::fabs(edt_bst(ipsi_gm_closed.vertices[n]));
            if (d_ipsi_gm_2_bst < 1.0f) {
                midLine.push_back(true);
                continue;
            }

            float d_ipsi_wm_2_vdc = std::fabs(edt_vdc(ipsi_wm_closed.vertices[n]));
            if (d_ipsi_wm_2_vdc < 1.0f) {
                midLine.push_back(true);
                continue;
            }

            midLine.push_back(false);            
        }

        SurfaceField midLineMask = ipsi_wm_closed.makeVertField("midLineMask", midLine);

        Surface ipsi_wm_open;
        removeVertices(&ipsi_wm_open,&ipsi_wm_closed, &midLineMask);
        ipsi_wm_open = surfMakeItSingleOpen(ipsi_wm_open);

        Surface ipsi_gm_open;
        removeVertices(&ipsi_gm_open,&ipsi_gm_closed, &midLineMask);
        ipsi_gm_open = surfMakeItSingleOpen(ipsi_gm_open);

        Surface gm_ribbon = surfGlueBoundaries(ipsi_wm_open, ipsi_gm_open);
        gm_ribbon = surfMakeItSingleClosed(gm_ribbon);
        gm_ribbon = surfRemesh(gm_ribbon,gm_ribbon.nv);
        gm_ribbon = surfMakeItSingleClosed(gm_ribbon);

        return gm_ribbon;

    };

    Surface l_gm_ribbon = genGMribbon(l_wm_closed,l_gm_closed,r_wm_closed); l_gm_ribbon.write(outputFolder + "/l_gm.vtk");
    Surface r_gm_ribbon = genGMribbon(r_wm_closed,r_gm_closed,l_wm_closed); r_gm_ribbon.write(outputFolder + "/r_gm.vtk");



    // Combined surfaces for tractography
    Surface wm = surfMerge(r_wm_closed, l_wm_closed);
    Surface gm = surfMerge(r_gm_ribbon, l_gm_ribbon);
    
    // SEED
    disp(MSG_INFO,"Generating a combined seed surface");
    Surface seed = surfMerge(wm, bst);
    seed         = surfMerge(seed, cer);
    seed.write(outputFolder + "/seed.vtk");
    disp(MSG_INFO,"seed saved.");

    // REQUIRE_END_INSIDE
    disp(MSG_INFO,"Generating a combined require_end_s inside surface");
    Surface req_end_inside  = surfMerge(gm, sgm);
    req_end_inside          = surfMerge(req_end_inside, cer);
    req_end_inside          = surfMerge(req_end_inside, bst);
    req_end_inside.write(outputFolder + "/require_end_inside.vtk");
    disp(MSG_INFO,"require_end_inside saved.");

    wm.write(outputFolder + "/wm.vtk");
    gm.write(outputFolder + "/gm.vtk");

    return;

    
    /*


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

    // Calculate distance images
    Image<float> edt_l_gm;
    edt_l_gm.createFromTemplate(ribbon,true);
    mapSurface2Image(&l_gm_closed, &edt_l_gm, 0, NULL, NULL, EDT);

    Image<float> edt_r_gm;
    edt_r_gm.createFromTemplate(ribbon,true);
    mapSurface2Image(&r_gm_closed, &edt_r_gm, 0, NULL, NULL, EDT);


    // Compute midline distance
    std::vector<int> l_midLine;
    l_midLine.reserve(l_wm_closed.nv);

    for (int n = 0; n < l_wm_closed.nv; n++) {

        bool isMidLine = true;

        float val = std::fabs(edt_r_gm(l_wm_closed.vertices[n]));
        isMidLine = (val < 1.0f);

        // for (const auto& p : nei) {
        //     float r[3];
        //     r[0] = p[0] + l_wm_closed.vertices[n][0];
        //     r[1] = p[1] + l_wm_closed.vertices[n][1];
        //     r[2] = p[2] + l_wm_closed.vertices[n][2];
        //     float val = gm_ribbon(r);
        //     if (val > EPS6) {
        //         isMidLine = false;
        //         break;
        //     }
        // }

        // if (isMidLine) {
        //     for (const auto& p : nei) {
        //         float r[3];
        //         r[0] = p[0] + l_gm_closed.vertices[n][0];
        //         r[1] = p[1] + l_gm_closed.vertices[n][1];
        //         r[2] = p[2] + l_gm_closed.vertices[n][2];
        //         float val = gm_ribbon(r);
        //         if (val > EPS6) {
        //             isMidLine = false;
        //             break;
        //         }
        //     }
        // }

        l_midLine.push_back(int(isMidLine));
    }

    SurfaceField l_midLineMask = l_wm_closed.makeVertField("midLineMask", l_midLine);

    Surface l_wm_open;
    removeVertices(&l_wm_open,&l_wm_closed, &l_midLineMask);
    l_wm_open.write(outputFolder + "/l_wm_open.vtk");

    Surface l_gm_open;
    removeVertices(&l_gm_open,&l_gm_closed, &l_midLineMask);
    // l_gm_open.write(outputFolder + "/l_gm_open.vtk");

    Surface l_gm_ribbon = surfGlueBoundaries(l_wm_open, l_gm_open);
    l_gm_ribbon = surfMakeItSingleClosed(l_gm_ribbon);
    l_gm_ribbon = surfRemesh(l_gm_ribbon,l_gm_ribbon.nv);
    l_gm_ribbon = surfMakeItSingleClosed(l_gm_ribbon);
    l_gm_ribbon.write(outputFolder + "/l_gm.vtk");
    


    //
    

    std::vector<int> r_midLine;
    r_midLine.reserve(r_wm_closed.nv);

    for (int n = 0; n < r_wm_closed.nv; n++) {

        bool isMidLine = true;

        float val = std::fabs(edt_l_gm(r_wm_closed.vertices[n]));
        isMidLine = (val < 1.0f);

        // for (const auto& p : nei) {
        //     float r[3];
        //     r[0] = p[0] + r_wm_closed.vertices[n][0];
        //     r[1] = p[1] + r_wm_closed.vertices[n][1];
        //     r[2] = p[2] + r_wm_closed.vertices[n][2];
        //     float val = gm_ribbon(r);
        //     if (val > EPS6) {
        //         isMidLine = false;
        //         break;
        //     }
        // }

        // if (isMidLine) {
        //     for (const auto& p : nei) {
        //         float r[3];
        //         r[0] = p[0] + r_gm_closed.vertices[n][0];
        //         r[1] = p[1] + r_gm_closed.vertices[n][1];
        //         r[2] = p[2] + r_gm_closed.vertices[n][2];
        //         float val = gm_ribbon(r);
        //         if (val > EPS6) {
        //             isMidLine = false;
        //             break;
        //         }
        //     }
        // }

        r_midLine.push_back(int(isMidLine));
    }

    SurfaceField r_midLineMask = r_wm_closed.makeVertField("midLineMask", r_midLine);

    Surface r_wm_open;
    removeVertices(&r_wm_open,&r_wm_closed, &r_midLineMask);
    r_wm_open.write(outputFolder + "/r_wm_open.vtk");

    Surface r_gm_open;
    removeVertices(&r_gm_open,&r_gm_closed, &r_midLineMask);
    // r_gm_open.write(outputFolder + "/r_gm_open.vtk");

    Surface r_gm_ribbon = surfGlueBoundaries(r_wm_open, r_gm_open);
    r_gm_ribbon = surfMakeItSingleClosed(r_gm_ribbon);
    r_gm_ribbon = surfRemesh(r_gm_ribbon,r_gm_ribbon.nv);
    r_gm_ribbon = surfMakeItSingleClosed(r_gm_ribbon);

    r_gm_ribbon.write(outputFolder + "/r_gm.vtk");


    */

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
