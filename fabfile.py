from fabric.api import *
from xml.etree import ElementTree as ET

env.hosts = ['proteomics2.ucsd.edu']
env.user = 'bpullman'
env.workflow_components = ['input.xml', 'binding.xml', 'flow.xml', 'result.xml', 'tool.xml']

def updateWorkflowComponent(workflow,component):
    local = '{}/{}'.format(workflow,component)
    xml = ET.parse(local)
    server = '/ccms/workflows/{}/{}'.format(workflow, component)
    put(local, server)

def updateWorkflow(workflow):
    for component in env.workflow_components:
        updateWorkflowComponent(workflow,component)

def updatePeptideStatisticsWorkflow():
    updateWorkflow('peptide_statistics_hpp')

def updatePeptideStatisticsScripts():
    put('peptide_statistics_hpp_tools/*', '/data/cluster/tools/peptide_statistics_hpp')

def updateAll():
    updatePeptideStatisticsWorkflow()
    updatePeptideStatisticsScripts()
