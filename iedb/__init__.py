# **************************************************************************
# *
# * Authors:	Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
# *			 	Carlos Oscar Sorzano (coss@cnb.csic.es)
# *			 	Martín Salinas Antón (martin.salinas@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# * All comments concerning this program package may be sent to the
# * e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""
This package contains protocols for creating and using ConPLex models for virtual screening
"""

# General imports
import os, subprocess

# Scipion em imports
from pwem import Config as emConfig

# Plugin imports
from pwchem import Plugin as pwchemPlugin
from pwchem.utils import insistentRun

from .bibtex import _bibtexStr
from .constants import *

# Pluging variables
_logo = 'dtu_logo.jpeg'

class Plugin(pwchemPlugin):
	"""
	"""

	@classmethod
	def _defineVariables(cls):
		#  BepiPred
		cls._defineVar(BEPIPRED_DIC['home'], cls.getDefaultDir(BEPIPRED_DIC))
		cls._defineVar(BEPIPRED_DIC['activation'], cls.getEnvActivationCommand(BEPIPRED_DIC))
		cls._defineVar(BEPIPRED_DIC['zip'], None)

		# MHC-I
		cls._defineVar(MHCI_DIC['home'], cls.getDefaultDir(MHCI_DIC))
		cls._defineVar(MHCI_DIC['tar'], None)

		# MHC-II
		cls._defineVar(MHCII_DIC['home'], cls.getDefaultDir(MHCII_DIC))
		cls._defineVar(MHCII_DIC['tar'], None)

		# Population Coverage
		cls._defineVar(COVE_DIC['home'], cls.getDefaultDir(COVE_DIC))
		cls._defineVar(COVE_DIC['tar'], None)

		cls._defineVar(ELLI_DIC['home'], cls.getDefaultDir(ELLI_DIC))
		cls._defineVar(ELLI_DIC['jar'], None)

		cls._defineVar(IMMU_DIC['home'], cls.getDefaultDir(IMMU_DIC))
		cls._defineVar(IMMU_DIC['tar'], None)
		# The vendored predict_immunogenicity.py is a Python 2 script; installed into its own
		# conda env (see _addImmunogenicityPackage) so users don't need to build one themselves.
		cls._defineVar(IMMU_DIC['activation'], cls.getEnvActivationCommand(IMMU_DIC))

		# IEDB B-cell protective/neutralizing epitope reference (manually
		# downloaded, see ProtIEDBCrossref/constants.py)
		cls._defineVar(IEDB_BCELL_REFERENCE_PATH, '')

	@classmethod
	def defineBinaries(cls, env):
		"""This function defines the binaries for each package."""
		#  B-CELL package
		if cls.checkVarPath(BEPIPRED_DIC, 'zip'):
			cls._addBepiPredPackage(env, zipPath=cls.getVar(BEPIPRED_DIC['zip']))
		elif cls.checkVarPath(BEPIPRED_DIC, 'home'):
			if not cls.checkCallEnv(BEPIPRED_DIC):
				cls._addBepiPredPackage(env, bepiHome=cls.getVar(BEPIPRED_DIC['home']))
			# else:
			# 	print('Environment activation command and HOME variables for BepiPred already found, installation no needed')

	  #  MHC-I package
		if cls.checkVarPath(MHCI_DIC, 'tar'):
			cls._addMHCIPackage(env, tarPath=cls.getVar(MHCI_DIC['tar']))
		elif cls.checkVarPath(MHCI_DIC, 'home'):
			cls._addMHCIPackage(env, mhciHome=cls.getVar(MHCI_DIC['home']))

		#  MHC-II package
		if cls.checkVarPath(MHCII_DIC, 'tar'):
			cls._addMHCIIPackage(env, tarPath=cls.getVar(MHCII_DIC['tar']))
		elif cls.checkVarPath(MHCII_DIC, 'home'):
			cls._addMHCIIPackage(env, mhciiHome=cls.getVar(MHCII_DIC['home']))

		#  Population Coverage package
		if cls.checkVarPath(COVE_DIC, 'tar'):
			cls._addCoveragePackage(env, tarPath=cls.getVar(COVE_DIC['tar']))
		elif cls.checkVarPath(COVE_DIC, 'home'):
			cls._addCoveragePackage(env, coveHome=cls.getVar(COVE_DIC['home']))

		if cls.checkVarPath(ELLI_DIC, 'jar'):
			cls._addElliProPackage(env, elliJar=cls.getVar(ELLI_DIC['jar']))

		#  Class I Immunogenicity package
		if cls.checkVarPath(IMMU_DIC, 'tar'):
			cls._addImmunogenicityPackage(env, tarPath=cls.getVar(IMMU_DIC['tar']))
		elif cls.checkVarPath(IMMU_DIC, 'home'):
			cls._addImmunogenicityPackage(env, immunoHome=cls.getVar(IMMU_DIC['home']))


	@classmethod
	def _addBepiPredPackage(cls, env, bepiHome=None, zipPath=None, default=True):
		""" This function provides the neccessary commands for installing BepiPred. """
		BEPIPRED_INSTALLED = '%s_installed' % BEPIPRED_DIC['name']

		installationCmd = ''
		if not bepiHome and zipPath:
			bepiHome = os.path.join(emConfig.EM_ROOT, cls.getEnvName(BEPIPRED_DIC))
			installationCmd += f'unzip -q {zipPath} -d {bepiHome} && ' \
												 f'mv {bepiHome}/BepiPred3_src/* {bepiHome} && rm -r {bepiHome}/BepiPred3_src && '

		installationCmd += f"cd {bepiHome} && sed -i 's/^torch==1.12.0/torch/g' requirements.txt && "
		installationCmd += f"conda create -y -n {cls.getEnvName(BEPIPRED_DIC)} python=3.8 && "
		installationCmd += f"{cls.getCondaActivationCmd()}conda activate {cls.getEnvName(BEPIPRED_DIC)} && pip install -r requirements.txt && "
		installationCmd += f"touch {BEPIPRED_INSTALLED}"

		env.addPackage(BEPIPRED_DIC['name'], version=BEPIPRED_DIC['version'],
									 commands=[(installationCmd, os.path.join(bepiHome, BEPIPRED_INSTALLED))], tar='void.tgz',
									 neededProgs=["conda"], default=default, buildDir=os.path.split(bepiHome)[-1])

	@classmethod
	def _addMHCIPackage(cls, env, mhciHome=None, tarPath=None, default=True):
		""" This function provides the neccessary commands for installing the MHC-I prediction program. """
		MHCI_INSTALLED = '%s_installed' % MHCI_DIC['name']
		emHome = os.path.join(emConfig.EM_ROOT, cls.getEnvName(MHCI_DIC))

		installationCmd = ''
		if not mhciHome and tarPath:
			mhciHome = emHome
			installationCmd += f'tar -xf {tarPath} -C {mhciHome} && ' \
												 f'mv {mhciHome}/mhc_i/* {mhciHome} && rm -r {mhciHome}/mhc_i && '

		if mhciHome != emHome:
			installationCmd += f"mv {mhciHome}/* {emHome} && rm -r {mhciHome} && "
		installationCmd += f"cd {emHome} && ./configure && "
		installationCmd += f"touch {MHCI_INSTALLED}"

		env.addPackage(MHCI_DIC['name'], version=MHCI_DIC['version'],
									commands=[(installationCmd, os.path.join(emHome, MHCI_INSTALLED))], tar='void.tgz',
									default=default, buildDir=os.path.split(mhciHome)[-1])

	@classmethod
	def _addMHCIIPackage(cls, env, mhciiHome=None, tarPath=None, default=True):
		""" This function provides the neccessary commands for installing the MHC-II prediction program. """
		MHCII_INSTALLED = '%s_installed' % MHCII_DIC['name']
		emHome = os.path.join(emConfig.EM_ROOT, cls.getEnvName(MHCII_DIC))

		installationCmd = ''
		if not mhciiHome and tarPath:
			mhciiHome = emHome
			installationCmd += f'tar -xf {tarPath} -C {mhciiHome} && ' \
												 f'mv {mhciiHome}/mhc_ii/* {mhciiHome} && rm -r {mhciiHome}/mhc_ii && '

		if mhciiHome != emHome:
			installationCmd += f"mv {mhciiHome}/* {emHome} && rm -r {mhciiHome} && "
		installationCmd += f"cd {emHome} && ./configure.py && "
		installationCmd += f"touch {MHCII_INSTALLED}"

		env.addPackage(MHCII_DIC['name'], version=MHCII_DIC['version'],
									 commands=[(installationCmd, os.path.join(emHome, MHCII_INSTALLED))], tar='void.tgz',
									 default=default, buildDir=os.path.split(mhciiHome)[-1])

	@classmethod
	def _addCoveragePackage(cls, env, coveHome=None, tarPath=None, default=True):
		""" This function provides the neccessary commands for installing the MHC coverage program. """
		COVE_INSTALLED = '%s_installed' % COVE_DIC['name']
		emHome = os.path.join(emConfig.EM_ROOT, cls.getEnvName(COVE_DIC))

		installationCmd = ''
		if not coveHome and tarPath:
			coveHome = emHome
			installationCmd += f'tar -xf {tarPath} -C {coveHome} && ' \
							   f'mv {coveHome}/population_coverage/* {coveHome} && rm -r {coveHome}/population_coverage && '

		if coveHome != emHome:
			installationCmd += f"mv {coveHome}/* {emHome} && rm -r {coveHome} && "
		installationCmd += f"cd {emHome} && ./configure && "
		installationCmd += f"touch {COVE_INSTALLED}"

		env.addPackage(COVE_DIC['name'], version=COVE_DIC['version'],
					   commands=[(installationCmd, os.path.join(emHome, COVE_INSTALLED))], tar='void.tgz',
					   default=default, buildDir=os.path.split(coveHome)[-1])

	@classmethod
	def _addElliProPackage(cls, env, elliJar, default=True):
		""" This function provides the neccessary commands for installing ElliPro. """
		ELLI_INSTALLED = '%s_installed' % ELLI_DIC['name']
		emHome = os.path.join(emConfig.EM_ROOT, cls.getEnvName(ELLI_DIC))

		installationCmd = f"mv {elliJar} {emHome} && cd {emHome} && touch {ELLI_INSTALLED}"
		env.addPackage(ELLI_DIC['name'], version=ELLI_DIC['version'],
									 commands=[(installationCmd, os.path.join(emHome, ELLI_INSTALLED))], tar='void.tgz',
									 default=default)

	@classmethod
	def _addImmunogenicityPackage(cls, env, immunoHome=None, tarPath=None, default=True):
		""" This function provides the neccessary commands for installing the Class I Immunogenicity program. """
		IMMUNO_INSTALLED = '%s_installed' % IMMU_DIC['name']
		emHome = os.path.join(emConfig.EM_ROOT, cls.getEnvName(IMMU_DIC))

		installationCmd = ''
		if not immunoHome and tarPath:
			immunoHome = emHome
			installationCmd += f'tar -xf {tarPath} -C {immunoHome} && ' \
												 f'mv {immunoHome}/immunogenicity/* {immunoHome} && rm -r {immunoHome}/immunogenicity && '

		if immunoHome != emHome:
			installationCmd += f"mv {immunoHome}/* {emHome} && rm -r {immunoHome} && "
		# The vendored predict_immunogenicity.py is a Python 2 script; create a dedicated conda
		# env for it here so users don't need any extra manual step.
		installationCmd += f"conda create -y -n {cls.getEnvName(IMMU_DIC)} python=2.7 && "
		installationCmd += f"touch {IMMUNO_INSTALLED}"

		env.addPackage(IMMU_DIC['name'], version=IMMU_DIC['version'],
									 commands=[(installationCmd, os.path.join(emHome, IMMUNO_INSTALLED))], tar='void.tgz',
									 neededProgs=["conda"], default=default, buildDir=os.path.split(immunoHome)[-1])


	@classmethod
	def validateBepiPredInstallation(cls):
		""" Check if the BepiPred installation alone is correct. Returning an empty list means that the
		installation is correct and there are not errors. If some errors are found, a list with the error
		messages will be returned."""
		mPaths = []
		if not cls.checkVarPath(BEPIPRED_DIC, 'home'):
			mPaths.append(f"Path of BepiPred home (folder like BepiPred3_src) does not exist.\n"
										f"You must either define it in the scipion.conf (as {BEPIPRED_DIC['home']} = <pathToBepiPred_src>) "
										f"or define the location of the raw dowloaded ZIP file (like bepipred-3.0b.src.zip) as "
										f"{BEPIPRED_DIC['zip']} = <pathToBepiPredZip>.\nAlternatively, you can move the home folder into "
										f"{emConfig.EM_ROOT} keeping the '{BEPIPRED_DIC['pattern']}' pattern.")
		elif not cls.checkCallEnv(BEPIPRED_DIC):
			mPaths.append(f"Activation of the BepiPred environment failed.\n")
		return mPaths

	@classmethod
	def validatePackageInstallation(cls, softDic, progFile):
		""" Check if a single manually-installed IEDB package (MHC-I, MHC-II, population coverage or
		immunogenicity) is correct. Returning an empty list means that the installation is correct and
		there are not errors. If some errors are found, a list with the error messages will be returned."""
		mPaths = []
		if not cls.checkVarPath(softDic, 'home') or \
				not os.path.exists(os.path.join(cls.getVar(softDic['home']), progFile)):
			mPaths.append(f"Path of {softDic['name']} home does not exist or is incomplete.\n"
										f"You must either define it in the scipion.conf (as {softDic['home']} = <pathTo{softDic['name']}Folder>) "
										f"or define the location of the raw downloaded tar file as "
										f"{softDic['tar']} = <pathTo{softDic['name']}Tar>.\nAlternatively, you can move the home folder into "
										f"{emConfig.EM_ROOT} keeping the '{softDic['pattern']}' pattern.")
		return mPaths

	@classmethod
	def validateMHCIInstallation(cls):
		return cls.validatePackageInstallation(MHCI_DIC, 'src/predict_binding.py')

	@classmethod
	def validateMHCIIInstallation(cls):
		return cls.validatePackageInstallation(MHCII_DIC, 'mhc_II_binding.py')

	@classmethod
	def validateCoverageInstallation(cls):
		return cls.validatePackageInstallation(COVE_DIC, 'calculate_population_coverage.py')

	@classmethod
	def validateImmunogenicityInstallation(cls):
		""" Check if the Immunogenicity installation, including its own Python 2 conda env,
		is correct. Returning an empty list means that the installation is correct and there
		are not errors. If some errors are found, a list with the error messages will be
		returned."""
		mPaths = cls.validatePackageInstallation(IMMU_DIC, 'predict_immunogenicity.py')
		if not mPaths and not cls.checkCallEnv(IMMU_DIC):
			mPaths.append("Activation of the Immunogenicity environment failed.\n")
		return mPaths

	@classmethod
	def validateElliProInstallation(cls):
		""" Check if the ElliPro installation alone is correct. Returning an empty list means that the
		installation is correct and there are not errors. If some errors are found, a list with the error
		messages will be returned."""
		mPaths = []
		if not cls.checkVarPath(ELLI_DIC, 'home') or \
				not os.path.exists(os.path.join(cls.getVar(ELLI_DIC['home']), 'ElliPro.jar')):
			mPaths.append(f"ElliPro.jar was not found in {ELLI_DIC['home']}.\n"
										f"You must provide the downloaded jar file location in the scipion.conf as "
										f"{ELLI_DIC['jar']} = <pathToElliProJar>.")
		return mPaths

	@classmethod
	def validateInstallation(cls):
		""" Check if the installation of the plugin as a whole is correct, across all six packages it
		provides. Returning an empty list means that the installation is correct and there are not
		errors. If some errors are found, a list with the error messages will be returned.
		This is the plugin-wide check used e.g. by the plugin manager; individual protocols override
		their own ``validateInstallation`` classmethod so that running one protocol only requires the
		package it actually needs, not every package this plugin bundles."""
		mPaths = cls.validateBepiPredInstallation() + cls.validateMHCIInstallation() + cls.validateMHCIIInstallation() + \
						 cls.validateCoverageInstallation() + cls.validateImmunogenicityInstallation() + cls.validateElliProInstallation()

		if len(mPaths) > 0:
			mPaths.append(NOINSTALL_WARNING)
		return mPaths

	@classmethod
	def getDefaultDir(cls, softDic, fn=""):
		emDir = emConfig.EM_ROOT
		pattern = softDic['pattern'].lower()
		for file in os.listdir(emDir):
			fileLower = file.lower()
			# The directory name must start with the pattern, and whatever follows it (if anything)
			# must be a separator or a digit, not a letter: a bare substring check would match
			# 'mhc_ii-...' against pattern 'mhc_i' (since 'mhc_i' is itself a substring of 'mhc_ii'),
			# silently resolving MHC-I's home to the MHC-II install. Allowing a digit right after the
			# pattern (no separator) keeps matching folder names like 'BepiPred3_src'.
			if fileLower.startswith(pattern):
				rest = fileLower[len(pattern):]
				if rest == '' or rest[0] in '-_' or rest[0].isdigit():
					foundDir = os.path.join(emDir, file, fn)
					return foundDir.rstrip('/')
		return os.path.join(emConfig.EM_ROOT, cls.getEnvName(softDic))

	@classmethod
	def checkVarPath(cls, softDic, var='home'):
		'''Check if a plugin variable exists and so do its path'''
		exists = False
		varValue = cls.getVar(softDic[var])
		if varValue and os.path.exists(varValue):
			exists = True
		return exists

	@classmethod
	def checkCallEnv(cls, packageDic):
		actCommand = cls.getVar(packageDic['activation'])
		try:
			if 'conda' in actCommand and not 'shell.bash hook' in actCommand:
				actCommand = f'{cls.getCondaActivationCmd()}{actCommand}'
			subprocess.check_output(actCommand, shell=True)
			envFine = True
		except subprocess.CalledProcessError as e:
			envFine = False
		return envFine

	@classmethod
	def getPluginHome(cls, path=""):
		import iedb
		fnDir = os.path.split(iedb.__file__)[0]
		return os.path.join(fnDir, path)

	# ---------------------------------- Protocol functions-----------------------
	@classmethod
	def runBepiPred(cls, protocol, args, cwd=None, popen=False):
		""" Run rdkit command from a given protocol. """
		bepiHome, bepiAct = cls.getVar(BEPIPRED_DIC["home"]), cls.getVar(BEPIPRED_DIC["activation"])
		fullProgram = f'{bepiAct} && python {os.path.join(bepiHome, "bepipred3_CLI.py")}'
		if not popen:
			protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)
		else:
			subprocess.check_call(f'{fullProgram} {args}', cwd=cwd, shell=True)

	@classmethod
	def runMHC_I(cls, protocol, args, cwd=None, popen=False):
		""" Run mhc-i command from a given protocol. """
		mhciHome = cls.getVar(MHCI_DIC["home"])
		fullProgram = f'{os.path.join(mhciHome, "src/predict_binding.py")}'
		if not popen:
			protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)
		else:
			subprocess.check_call(f'{fullProgram} {args}', cwd=cwd, shell=True)

	@classmethod
	def runMHC_II(cls, protocol, args, cwd=None, popen=False):
		""" Run mhc-ii command from a given protocol. """
		mhciiHome = cls.getVar(MHCII_DIC["home"])
		fullProgram = f'python {os.path.join(mhciiHome, "mhc_II_binding.py")}'
		if not popen:
			protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)
		else:
			subprocess.check_call(f'{fullProgram} {args}', cwd=cwd, shell=True)

	@classmethod
	def runPopulationCoverage(cls, args, protocol=None, cwd=None, popen=False):
		""" Run population coverage command from a given protocol. """
		coveHome = cls.getVar(COVE_DIC["home"])
		fullProgram = f'python {os.path.join(coveHome, "calculate_population_coverage.py")} '
		print('Running coverage: ', fullProgram + args)
		if not popen and protocol:
			insistentRun(protocol, fullProgram, args, kwargs={'cwd': cwd})
			# protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)
		else:
			subprocess.check_call(f'{fullProgram} {args}', cwd=cwd, shell=True)

	@classmethod
	def runElliPro(cls, protocol, args, cwd=None, popen=False):
		""" Run ElliPro command from a given protocol. """
		elliHome = cls.getVar(ELLI_DIC["home"])
		fullProgram = f'java -jar {os.path.join(elliHome, "ElliPro.jar")}'
		if not popen:
			protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)
		else:
			subprocess.check_call(f'{fullProgram} {args}', cwd=cwd, shell=True)

	@classmethod
	def runImmunogenicity(cls, protocol, args, cwd=None, popen=False):
		""" Run immunogenicity command from a given protocol. """
		immuHome, immuAct = cls.getVar(IMMU_DIC["home"]), cls.getVar(IMMU_DIC["activation"])
		fullProgram = f'{immuAct} && python {os.path.join(immuHome, "predict_immunogenicity.py")}'
		if not popen:
			protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)
		else:
			subprocess.check_call(f'{fullProgram} {args}', cwd=cwd, shell=True)

	# ---------------------------------- Utils functions-----------------------

