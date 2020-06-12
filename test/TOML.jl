using TOML

@testset "TOML.jl" begin
	params = TOML.parsefile("TOML.toml")
	
	# test getTable function
	@test MPISimulations.getTable(params,"Test.Table") == Dict{AbstractString,Any}()
	@test_throws MissingException MPISimulations.getTable(params,"missingTable")
	@test ismissing(MPISimulations.getTable(params,"missingTable",tableOptional=true))

	#test getValue function
	@test MPISimulations.getValue(params,"Test","key") == 42
	@test_throws MissingException  MPISimulations.getValue(params,"Test","missingKey")
	@test_throws MissingException  MPISimulations.getValue(params,"missingTable","missingKey")
	@test ismissing(MPISimulations.getValue(params,"Test","missingKey",valueOptional=true))
	@test ismissing(MPISimulations.getValue(params,"missingTable","missingKey",valueOptional=true))
	
	# test getType function
	@test MPISimulations.getType(params,"Test") == ComplexF64
	@test ismissing(MPISimulations.getType(params,"Test.Table"))
	@info "" MPISimulations.getType(params,"Test.Table")

	# test initialize function
	@test MPISimulations.initialize(Complex,params,Float64) == ComplexF64(1.0,0.5)
	@test MPISimulations.initialize(Complex,params,Float32) == ComplexF32(1.0,0.5)

	# test getTableName function
	@test MPISimulations.getTableName(Complex) == "Number.Complex"
end
