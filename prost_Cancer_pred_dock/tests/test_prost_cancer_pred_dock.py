from prost_cancer_pred_dock import __version__
from prost_cancer_pred_dock import predict_pIC50

def test_version():
    assert __version__ == '0.0.1'

def test_predict_pIC50():
    results = predict_pIC50("Oc1ccc2c(c1)S[C@H](c1ccco1)[C@H](c1ccc(OCCN3CCCCC3)cc1)O2")
    assert results == '7.62'