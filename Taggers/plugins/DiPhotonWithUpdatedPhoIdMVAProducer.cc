
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "flashgg/MicroAOD/interface/PhotonIdUtils.h"
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/DiPhotonPhoIdMVAInputVars.h"

#include "TFile.h"
#include "TGraph.h"

namespace flashgg {

    class DiPhotonWithUpdatedPhoIdMVAProducer : public edm::EDProducer
    {
    public:
        DiPhotonWithUpdatedPhoIdMVAProducer( const edm::ParameterSet & );
        void produce( edm::Event &, const edm::EventSetup & ) override;

    private:
        edm::EDGetTokenT<edm::View<flashgg::DiPhotonCandidate> > token_;
        edm::EDGetTokenT<double> rhoToken_;
        PhotonIdUtils phoTools_;
        edm::FileInPath phoIdMVAweightfileEB_, phoIdMVAweightfileEE_, correctionFile_;
        bool correctInputs_;
        bool debug_;
        //        std::vector<TGraph*> corrections_;
        std::vector<std::unique_ptr<TGraph> > corrections_;
    };

    DiPhotonWithUpdatedPhoIdMVAProducer::DiPhotonWithUpdatedPhoIdMVAProducer( const edm::ParameterSet &ps ) :
        token_(consumes<edm::View<flashgg::DiPhotonCandidate> >(ps.getParameter<edm::InputTag>("src"))),
        rhoToken_( consumes<double>( ps.getParameter<edm::InputTag>( "rhoFixedGridCollection" ) ) ),
        debug_( ps.getParameter<bool>( "Debug" ) )
    {
        phoIdMVAweightfileEB_ = ps.getParameter<edm::FileInPath>( "photonIdMVAweightfile_EB" );
        phoIdMVAweightfileEE_ = ps.getParameter<edm::FileInPath>( "photonIdMVAweightfile_EE" );
        phoTools_.setupMVA( phoIdMVAweightfileEB_.fullPath(), phoIdMVAweightfileEE_.fullPath() );

        correctInputs_ = ps.existsAs<edm::FileInPath>("correctionFile") ? true: false;
        if (correctInputs_) {
            correctionFile_ = ps.getParameter<edm::FileInPath>( "correctionFile" );
            TFile* f = TFile::Open(correctionFile_.fullPath().c_str());
            corrections_.emplace_back((TGraph*)((TGraph*) f->Get("transffull5x5R9EB"))->Clone() );
            corrections_.emplace_back((TGraph*)((TGraph*) f->Get("transfEtaWidthEB"))->Clone() );
            corrections_.emplace_back((TGraph*)((TGraph*) f->Get("transfS4EB"))->Clone() );
            corrections_.emplace_back((TGraph*)((TGraph*) f->Get("transffull5x5R9EE"))->Clone() );
            corrections_.emplace_back((TGraph*)((TGraph*) f->Get("transfEtaWidthEE"))->Clone() );
            corrections_.emplace_back((TGraph*)((TGraph*) f->Get("transfS4EE"))->Clone() );
            f->Close();
        }

        produces<std::vector<flashgg::DiPhotonCandidate> >();

        // map of diphoton candidate to values of input variables for Photon ID MVA
        produces<flashgg::DiPhotonPhoIdMVAInputVarsAssociation>();
    }

    void DiPhotonWithUpdatedPhoIdMVAProducer::produce( edm::Event &evt, const edm::EventSetup & )
    {
        edm::Handle<edm::View<flashgg::DiPhotonCandidate> > objects;
        evt.getByToken( token_, objects );

        edm::Handle<double> rhoHandle;
        evt.getByToken( rhoToken_, rhoHandle );
        const double rhoFixedGrd = *( rhoHandle.product() );

        auto_ptr<std::vector<flashgg::DiPhotonCandidate> > out_obj( new std::vector<flashgg::DiPhotonCandidate>() );
        // allocate objects for all photons now already
        std::vector<DiPhotonPhoIdMVAInputVars> phoIdMVAInputVars(objects->size());

        for (const auto & obj : *objects) {
            flashgg::DiPhotonCandidate *new_obj = obj.clone();
            new_obj->makePhotonsPersistent();
            double leadCorrectedEtaWidth = 0., subLeadCorrectedEtaWidth = 0.;
            if (not evt.isRealData() and correctInputs_) { 
                if (new_obj->getLeadingPhoton().isEB()) {
                    if (this->debug_) {
                        std::cout << new_obj->getLeadingPhoton().full5x5_r9() << std::endl;
                        std::cout << new_obj->getLeadingPhoton().r9() << std::endl;
                    }
                    reco::Photon::ShowerShape newShowerShapes = new_obj->getLeadingPhoton().full5x5_showerShapeVariables();
                    newShowerShapes.e3x3 = corrections_[0]->Eval(new_obj->getLeadingPhoton().full5x5_r9())*new_obj->getLeadingPhoton().superCluster()->rawEnergy();
                    new_obj->getLeadingPhoton().full5x5_setShowerShapeVariables(newShowerShapes);
                    leadCorrectedEtaWidth = corrections_[1]->Eval(new_obj->getLeadingPhoton().superCluster()->etaWidth());
                    new_obj->getLeadingPhoton().getSuperCluster()->setEtaWidth(leadCorrectedEtaWidth);
                    new_obj->getLeadingPhoton().setS4(corrections_[2]->Eval(new_obj->getLeadingPhoton().s4()));

                    if (this->debug_) {
                        std::cout << new_obj->getLeadingPhoton().full5x5_r9() << std::endl;
                        std::cout << new_obj->getLeadingPhoton().r9() << std::endl;
                    }
                }
                
                if (new_obj->getSubLeadingPhoton().isEB()) {
                    reco::Photon::ShowerShape newShowerShapes = new_obj->getSubLeadingPhoton().full5x5_showerShapeVariables();
                    newShowerShapes.e3x3 = corrections_[0]->Eval(new_obj->getSubLeadingPhoton().full5x5_r9())*new_obj->getSubLeadingPhoton().superCluster()->rawEnergy();
                    new_obj->getSubLeadingPhoton().full5x5_setShowerShapeVariables(newShowerShapes);
                    subLeadCorrectedEtaWidth = corrections_[1]->Eval(new_obj->getSubLeadingPhoton().superCluster()->etaWidth());
                    new_obj->getSubLeadingPhoton().getSuperCluster()->setEtaWidth(subLeadCorrectedEtaWidth);
                    new_obj->getSubLeadingPhoton().setS4(corrections_[2]->Eval(new_obj->getSubLeadingPhoton().s4()));
                }

                if (new_obj->getLeadingPhoton().isEE()) {
                    if (this->debug_) {
                        std::cout << new_obj->getLeadingPhoton().full5x5_r9() << std::endl;
                        std::cout << new_obj->getLeadingPhoton().r9() << std::endl;
                    }
                    reco::Photon::ShowerShape newShowerShapes = new_obj->getLeadingPhoton().full5x5_showerShapeVariables();
                    newShowerShapes.e3x3 = corrections_[3]->Eval(new_obj->getLeadingPhoton().full5x5_r9())*new_obj->getLeadingPhoton().superCluster()->rawEnergy();
                    new_obj->getLeadingPhoton().full5x5_setShowerShapeVariables(newShowerShapes);
                    leadCorrectedEtaWidth = corrections_[4]->Eval(new_obj->getLeadingPhoton().superCluster()->etaWidth());
                    new_obj->getLeadingPhoton().getSuperCluster()->setEtaWidth(leadCorrectedEtaWidth);
                    new_obj->getLeadingPhoton().setS4(corrections_[5]->Eval(new_obj->getLeadingPhoton().s4()));

                    if (this->debug_) {
                        std::cout << new_obj->getLeadingPhoton().full5x5_r9() << std::endl;
                        std::cout << new_obj->getLeadingPhoton().r9() << std::endl;
                    }
                }
                
                if (new_obj->getSubLeadingPhoton().isEE()) {
                    reco::Photon::ShowerShape newShowerShapes = new_obj->getSubLeadingPhoton().full5x5_showerShapeVariables();
                    newShowerShapes.e3x3 = corrections_[3]->Eval(new_obj->getSubLeadingPhoton().full5x5_r9())*new_obj->getSubLeadingPhoton().superCluster()->rawEnergy();
                    new_obj->getSubLeadingPhoton().full5x5_setShowerShapeVariables(newShowerShapes);
                    subLeadCorrectedEtaWidth = corrections_[4]->Eval(new_obj->getSubLeadingPhoton().superCluster()->etaWidth());
                    new_obj->getSubLeadingPhoton().getSuperCluster()->setEtaWidth(subLeadCorrectedEtaWidth);
                    new_obj->getSubLeadingPhoton().setS4(corrections_[5]->Eval(new_obj->getSubLeadingPhoton().s4()));
                }
            }

            if (this->debug_) {
                std::cout << " Input DiPhoton lead (sublead) MVA: " << obj.leadPhotonId() << " " << obj.subLeadPhotonId() << std::endl;
            }

            DiPhotonPhoIdMVAInputVars &thisPhoIdMVAInputVars = phoIdMVAInputVars[out_obj->size()];
            
            // calculate photon ID mva output, keep track of the values of the input variables
            float newleadmva = phoTools_.computeMVAWrtVtx( new_obj->getLeadingPhoton(), new_obj->vtx(), rhoFixedGrd, leadCorrectedEtaWidth, &thisPhoIdMVAInputVars.getInputsLeading());
            new_obj->getLeadingPhoton().setPhoIdMvaWrtVtx( new_obj->vtx(), newleadmva);
            float newsubleadmva = phoTools_.computeMVAWrtVtx( new_obj->getSubLeadingPhoton(), new_obj->vtx(), rhoFixedGrd, subLeadCorrectedEtaWidth, &thisPhoIdMVAInputVars.getInputsSubLeading());
            new_obj->getSubLeadingPhoton().setPhoIdMvaWrtVtx( new_obj->vtx(), newsubleadmva);
            if (this->debug_) {
                std::cout << " Output DiPhoton lead (sublead) MVA: " << new_obj->leadPhotonId() << " " << new_obj->subLeadPhotonId() << std::endl;
            }
            out_obj->push_back(*new_obj);
            delete new_obj;
        }

        // evt.put(..) returns an OrphanHandle (according to https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideEDMRef )
        // which does not have provenance information
        edm::OrphanHandle<std::vector<flashgg::DiPhotonCandidate> > diphotonCollHandle = evt.put(out_obj);

        // get edm::Ptrs back from the newly produced collection
        // (TODO: it would probably be cleaner to split this module into two: one which calculates
        //        a product with the photon id input variables (per diphoton) and one which
        //        then produces the updated diphotons with the new photon id output)
        auto_ptr<flashgg::DiPhotonPhoIdMVAInputVarsAssociation> inputMap(new flashgg::DiPhotonPhoIdMVAInputVarsAssociation());

        for (unsigned i = 0; i < diphotonCollHandle->size(); ++i)
        {
            inputMap->insert(make_pair(edm::Ptr<flashgg::DiPhotonCandidate>(diphotonCollHandle, i),
                                      phoIdMVAInputVars[i]));
        }

        evt.put(inputMap);
        
    }
}

typedef flashgg::DiPhotonWithUpdatedPhoIdMVAProducer FlashggDiPhotonWithUpdatedPhoIdMVAProducer;
DEFINE_FWK_MODULE( FlashggDiPhotonWithUpdatedPhoIdMVAProducer );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
