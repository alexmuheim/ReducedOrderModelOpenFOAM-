Reduced Order Model implementation for OpenFOAMv2312

To have more insights on how the code is used, refer to the READ_ME in the testCase folder.
Content:

pimpleParametric --> code to perform Full Order Model simulations for multiple viscosities. It automatically stores the fields with different names in the same OpenFOAM case folder.
					 
PODpost			 --> code to perform Proper Orthogonal Decomposition on pressure and velocity fields. 

offlineComputation --> performs the offline coefficient computation. It computes the tensors that are stored in the offline stage with different methods.  

galerkinProjection --> online computation for the Reduced Order Model. Creates residuals for the different stabilization methods and solves the ODE system with the implemented Newton solver. 
					   ResidualConstruction2eqs.H creates a DSL (Domain Specific Language) to automatically write the residuals in math form. 
 
errorComputation --> performs L2 relative error calculation as postProcessing. 
fieldDifference  --> computes the difference of 2 fields and stores them as a new field.

dictionaries --> all the relevant dictionaries needed to make the code work. They need to be placed in system apart for parameterDict which needs to be in constant. Look in testCase for more details.

testCase --> example case and usage instruction for the ROM. 
