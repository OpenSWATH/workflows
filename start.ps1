# workaround to keep window open
param($Work)

# restart PowerShell with -noexit, the same script, and 1
if (!$Work) {
    powershell -noexit -file $MyInvocation.MyCommand.Path 1
    return
}

# update docker with latest version
docker pull openswath/develop:latest

# start docker container mapping to current directory and remove image after execution
docker run --name osw --rm -v ${PWD}:/data -i -t openswath/develop:latest

