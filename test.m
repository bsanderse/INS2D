%TEST Run some test cases. Check that errors are bounded.

for testcase = {'LDC' 'BFS' 'shear_layer'}
    [errorV, errorp] = run_testsuite(testcase);
    assert(errorV < 1e-3);
end
