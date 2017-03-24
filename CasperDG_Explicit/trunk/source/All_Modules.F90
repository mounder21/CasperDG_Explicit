Module All_Modules
    use my_kinddefs
    use LU_module
    use inputs_module           !(dataModules)
    use norm_module
    use read_mesh_module           
    use constructGridData_module
    use get_h_module
    use initializeBasis_module
    use assembleLocMassMatrices_module
    use get_delta_t_module
    use allocate_module
    use residual_Jac_module
    use residual_module
    use computeJ_detJ_module 
    use initializeSolution_module
    use computeNormal_module
    use vtu_output_module
    use evolveSolution_module
    use solutionNorm_module
    use steadySolve_module
    use testing_module
end module All_Modules