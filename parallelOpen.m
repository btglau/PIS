% a script to correctly open the parallel pool
if ~strcmp(getenv('SLURM_JOB_ID'),'')
    %If running on cluster
    
    % delete any pre-existing pools
    delete(gcp('nocreate'))
    
    %Change location of job storage location so multiple jobs can run
    %simultaneously. If the folders don't get deleted, just remove the
    %jobs folders for completed or cancelled jobs.
        fprintf('\nInitializing local scheduler data...\n')
    currClust = parcluster('local');
        fprintf('DataLocation is %s\n',currClust.JobStorageLocation);
    local_clust_data = fullfile(getenv('SCRATCH'), getenv('SLURM_JOB_ID'));
    mkdir(local_clust_data);
    currClust.JobStorageLocation = local_clust_data;
        fprintf('New data loc is: %s\n',currClust.JobStorageLocation);

    %Before creating pool, check the number of allocated processors first
    allocatedProc = str2double(getenv('SLURM_CPUS_ON_NODE'));
    %For some reason, rusty/flatiron does not want a pool with 40 workers
    %(see https://www.mathworks.com/support/bugreports/details/2029922)
    if allocatedProc == 40;allocatedProc = allocatedProc - 2;end
        fprintf('Looks like %g processors allocated...\n',allocatedProc)
    
    %Open parallel processing pool and make note of the temp directory
    randWaitTime = 1+30*rand();%Multiple parallel pools starting on the same machine and time can cause problems...
        fprintf('Waiting %2.3f seconds to prevent parpool opening errors...\n',randWaitTime)
    pause(randWaitTime); 
    
    % parpool - do NOT time out the pool in the case of long calculations
    %parpool(currClust,min([allocatedProc maxCores]),'IdleTimeout',Inf);
    parpool(currClust,allocatedProc);
end