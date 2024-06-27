# coding: latin-1
# **************************************************************************
# *
# * Authors:  Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
# *
# * Biocomputing Unit, CNB-CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

# Common constants
DEFAULT_VERSION = '3.0'

# Package dictionaries
BEPIPRED_DIC =  {'name': 'BepiPred',    'version': '3.0', 'pattern': 'bepipred',
                 'home': 'BEPIPRED_HOME', 'zip': 'BEPIPRED_ZIP', 'activation': 'BEPIPRED_ACTIVATION_CMD'}

MHCI_DIC = {'name': 'mhc_i',    'version': '3.1.5', 'pattern': 'mhc_i',
            'home': 'MHC_I_HOME', 'tar': 'MHC_I_TAR'}

MHCII_DIC = {'name': 'mhc_ii',    'version': '3.1.11', 'pattern': 'mhc_ii',
             'home': 'MHC_II_HOME', 'tar': 'MHC_II_TAR'}

COVE_DIC = {'name': 'population_coverage',    'version': '3.0.2', 'pattern': 'population_coverage',
            'home': 'COVERAGE_HOME', 'tar': 'COVERAGE_TAR'}

ELLI_DIC = {'name': 'ellipro',    'version': '1.0', 'pattern': 'ellipro',
            'home': 'ELLIPRO_HOME', 'jar': 'ELLIPRO_JAR'}

IMMU_DIC = {'name': 'immunogenicity',    'version': '1.1', 'pattern': 'immunogenicity',
            'home': 'IMMUNO_HOME', 'tar': 'IMMUNO_TAR'}

READ_URL = 'https://github.com/scipion-chem/scipion-chem-bepipred'

NOINSTALL_WARNING = f'Installation could not be completed because BepiPred download or installation has not been found.\n' \
                    f'Please check the scipion-chem-iedb README file to see more details about how to proceed with ' \
                    f'the installation. You can find this instruction in {READ_URL}'


AB27 = ["HLA-A*01:01", "HLA-A*02:01", "HLA-A*02:03", "HLA-A*02:06", "HLA-A*03:01", "HLA-A*11:01", "HLA-A*23:01", "HLA-A*24:02", "HLA-A*26:01", "HLA-A*30:01",
        "HLA-A*30:02", "HLA-A*31:01", "HLA-A*32:01", "HLA-A*33:01", "HLA-A*68:01", "HLA-A*68:02", "HLA-B*07:02", "HLA-B*08:01", "HLA-B*15:01", "HLA-B*35:01",
        "HLA-B*40:01", "HLA-B*44:02", "HLA-B*44:03", "HLA-B*51:01", "HLA-B*53:01", "HLA-B*57:01", "HLA-B*58:01"]
MHCI_FREQ = ["HLA-A*01:01", "HLA-A*02:01", "HLA-A*02:06", "HLA-A*03:01", "HLA-A*11:01", "HLA-A*23:01", "HLA-A*24:02", "HLA-A*25:01", "HLA-A*26:01",
             "HLA-A*29:02", "HLA-A*30:01", "HLA-A*30:02", "HLA-A*31:01", "HLA-A*32:01", "HLA-A*33:01", "HLA-A*33:03", "HLA-A*68:01", "HLA-A*68:02",
             "HLA-A*74:01", "HLA-B*07:02", "HLA-B*08:01", "HLA-B*13:01", "HLA-B*13:02", "HLA-B*14:02", "HLA-B*15:01", "HLA-B*15:02", "HLA-B*15:25",
             "HLA-B*18:01", "HLA-B*27:02", "HLA-B*27:05", "HLA-B*35:01", "HLA-B*35:03", "HLA-B*37:01", "HLA-B*38:01", "HLA-B*39:01", "HLA-B*40:01",
             "HLA-B*40:02", "HLA-B*44:02", "HLA-B*44:03", "HLA-B*46:01", "HLA-B*48:01", "HLA-B*49:01", "HLA-B*50:01", "HLA-B*51:01", "HLA-B*52:01",
             "HLA-B*53:01", "HLA-B*55:01", "HLA-B*56:01", "HLA-B*57:01", "HLA-B*58:01", "HLA-B*58:02", "HLA-C*01:02", "HLA-C*02:02", "HLA-C*02:09",
             "HLA-C*03:02", "HLA-C*03:03", "HLA-C*03:04", "HLA-C*04:01", "HLA-C*05:01", "HLA-C*06:02", "HLA-C*07:01", "HLA-C*07:02", "HLA-C*07:04",
             "HLA-C*08:01", "HLA-C*08:02", "HLA-C*12:02", "HLA-C*12:03", "HLA-C*14:02", "HLA-C*15:02", "HLA-C*16:01", "HLA-C*17:01", "HLA-E*01:01",
             "HLA-E*01:03", "HLA-G*01:01", "HLA-G*01:02", "HLA-G*01:03", "HLA-G*01:04", "HLA-G*01:06"]
MHCI_REP = ["HLA-A*01:01", "HLA-A*02:01", "HLA-A*03:01", "HLA-A*24:02", "HLA-A*26:01", "HLA-B*07:02", "HLA-B*08:01", "HLA-B*27:05", "HLA-B*39:01",
            "HLA-B*40:01", "HLA-B*58:01", "HLA-B*15:01"]
MHCI_alleles_dic = {'Most frequent A, B': AB27, 'Frequent (>1%)': MHCI_FREQ, "Representative HLA supertypes": MHCI_REP}

DR7 = ['DRB1*03:01', 'DRB1*07:01', 'DRB1*15:01', 'DRB3*01:01', 'DRB3*02:02', 'DRB4*01:01', 'DRB5*01:01']
MHCII_FREQ = ['DRB1*01:01', 'DRB1*03:01', 'DRB1*04:01', 'DRB1*04:05', 'DRB1*07:01', 'DRB1*08:02', 'DRB1*09:01',
              'DRB1*11:01', 'DRB1*12:01', 'DRB1*13:02', 'DRB1*15:01', 'DRB3*01:01', 'DRB3*02:02', 'DRB4*01:01',
              'DRB5*01:01', 'DPA1*01/DPB1*04:01', 'DPA1*01:03/DPB1*02:01', 'DPA1*02:01/DPB1*01:01',
              'DPA1*02:01/DPB1*05:01', 'DPA1*03:01/DPB1*04:02', 'DQA1*01:01/DQB1*05:01', 'DQA1*01:02/DQB1*06:02',
              'DQA1*03:01/DQB1*03:02', 'DQA1*04:01/DQB1*04:02', 'DQA1*05:01/DQB1*02:01', 'DQA1*05:01/DQB1*03:01']
MHCII_alleles_dic = {'7-allele method': DR7, 'Most frequent 26': MHCII_FREQ}

MOUSE_MHCII_ALLELES = ['H2-IAb', 'H2-IAd', 'H2-IEd']


POP_DIC = {'Area': {'East Asia': {'Japan': ['Japan Oriental'], 'Korea; South': ['Korea; South Oriental'], 'Mongolia': ['Mongolia Oriental']},
                    'Northeast Asia': {'China': ['China Oriental'], 'Hong Kong': ['Hong Kong Oriental']},
                    'South Asia': {'Sri Lanka': ['Sri Lanka Asian'], 'Pakistan': ['Pakistan Asian', 'Pakistan Mixed'], 'India': ['India Asian']},
                    'Southeast Asia': {'Indonesia': ['Indonesia Austronesian'], 'Singapore': ['Singapore Austronesian', 'Singapore Oriental'], 'Philippines': ['Philippines Austronesian'], 'Borneo': ['Borneo Austronesian'], 'Malaysia': ['Malaysia Austronesian', 'Malaysia Oriental'], 'Vietnam': ['Vietnam Oriental'], 'Thailand': ['Thailand Oriental'], 'Taiwan': ['Taiwan Oriental']},
                    'Southwest Asia': {'Israel': ['Israel Arab', 'Israel Jew'], 'Iran': ['Iran Kurd', 'Iran Persian'], 'Oman': ['Oman Arab'], 'Lebanon': ['Lebanon Arab', 'Lebanon Mixed'], 'Saudi Arabia': ['Saudi Arabia Arab'], 'Jordan': ['Jordan Arab'], 'United Arab Emirates': ['United Arab Emirates Arab']},
                    'Europe': {'Turkey': ['Turkey Caucasoid'], 'Italy': ['Italy Caucasoid'], 'Czech Republic': ['Czech Republic Caucasoid', 'Czech Republic Other'], 'Scotland': ['Scotland Caucasoid'], 'France': ['France Caucasoid'], 'Georgia': ['Georgia Caucasoid', 'Georgia Kurd'], 'Slovakia': ['Slovakia Caucasoid'], 'Ireland South': ['Ireland South Caucasoid'], 'Wales': ['Wales Caucasoid'], 'Norway': ['Norway Caucasoid'], 'Belarus': ['Belarus Caucasoid'], 'Slovenia': ['Slovenia Caucasoid'], 'Germany': ['Germany Caucasoid'], 'Belgium': ['Belgium Caucasoid'], 'Spain': ['Spain Caucasoid', 'Spain Jew', 'Spain Other'], 'Ukraine': ['Ukraine Caucasoid'], 'Ireland Northern': ['Ireland Northern Caucasoid'], 'Netherlands': ['Netherlands Caucasoid'], 'Denmark': ['Denmark Caucasoid'], 'Poland': ['Poland Caucasoid'], 'Finland': ['Finland Caucasoid'], 'Macedonia': ['Macedonia Caucasoid'], 'Sweden': ['Sweden Caucasoid'], 'Croatia': ['Croatia Caucasoid'], 'Switzerland': ['Switzerland Caucasoid'], 'Russia': ['Russia Caucasoid', 'Russia Mixed', 'Russia Other', 'Russia Siberian'], 'Bulgaria': ['Bulgaria Caucasoid', 'Bulgaria Other'], 'Romania': ['Romania Caucasoid'], 'England': ['England Caucasoid', 'England Jew', 'England Mixed'], 'Portugal': ['Portugal Caucasoid'], 'Serbia': ['Serbia Caucasoid'], 'United Kingdom': ['United Kingdom Caucasoid'], 'Austria': ['Austria Caucasoid'], 'Greece': ['Greece Caucasoid']},
                    'East Africa': {'Kenya': ['Kenya Black'], 'Zambia': ['Zambia Black'], 'Zimbabwe': ['Zimbabwe Black'], 'Uganda': ['Uganda Black']},
                    'West Africa': {'Ghana': ['Ghana Black'], 'Liberia': ['Liberia Black'], 'Ivory Coast': ['Ivory Coast Black'], 'Gambia': ['Gambia Black'], 'Cape Verde': ['Cape Verde Black'], 'Senegal': ['Senegal Black'], 'Guinea-Bissau': ['Guinea-Bissau Black'], 'Nigeria': ['Nigeria Black'], 'Burkina Faso': ['Burkina Faso Black']},
                    'Central Africa': {'Sao Tome and Principe': ['Sao Tome and Principe Black'], 'Rwanda': ['Rwanda Black'], 'Equatorial Guinea': ['Equatorial Guinea Black'], 'Central African Republic': ['Central African Republic Black'], 'Congo': ['Congo Black'], 'Gabon': ['Gabon Black'], 'Cameroon': ['Cameroon Black']},
                    'North Africa': {'Sudan': ['Sudan Arab', 'Sudan Black', 'Sudan Mixed'], 'Algeria': ['Algeria Arab'], 'Morocco': ['Morocco Arab', 'Morocco Caucasoid'], 'Tunisia': ['Tunisia Arab', 'Tunisia Berber'], 'Ethiopia': ['Ethiopia Black'], 'Mali': ['Mali Black']},
                    'South Africa': {'South Africa': ['South Africa Black', 'South Africa Other']},
                    'West Indies': {'Jamaica': ['Jamaica Black'], 'Martinique': ['Martinique Black'], 'Cuba': ['Cuba Caucasoid', 'Cuba Mixed', 'Cuba Mulatto'], 'Trinidad and Tobago': ['Trinidad and Tobago Asian']},
                    'North America': {'Canada': ['Canada Amerindian'], 'United States': ['United States Amerindian', 'United States Asian', 'United States Austronesian', 'United States Black', 'United States Caucasoid', 'United States Hispanic', 'United States Mestizo', 'United States Polynesian'], 'Mexico': ['Mexico Amerindian', 'Mexico Mestizo']},
                    'Central America': {'Guatemala': ['Guatemala Amerindian'], 'Costa Rica': ['Costa Rica Mestizo']},
                    'South America': {'Brazil': ['Brazil Amerindian', 'Brazil Caucasoid', 'Brazil Mixed', 'Brazil Mulatto', 'Brazil Other'], 'Chile': ['Chile Amerindian', 'Chile Hispanic', 'Chile Mixed'], 'Peru': ['Peru Amerindian', 'Peru Mestizo'], 'Colombia': ['Colombia Amerindian', 'Colombia Black', 'Colombia Mestizo'], 'Paraguay': ['Paraguay Amerindian'], 'Argentina': ['Argentina Amerindian', 'Argentina Caucasoid'], 'Bolivia': ['Bolivia Amerindian'], 'Venezuela': ['Venezuela Amerindian', 'Venezuela Caucasoid', 'Venezuela Mestizo', 'Venezuela Mixed'], 'Ecuador': ['Ecuador Amerindian', 'Ecuador Black']},
                    'Oceania': {'Tonga': ['Tonga Polynesian'], 'Australia': ['Australia Australian Aborigines', 'Australia Caucasoid'], 'Samoa': ['Samoa Polynesian'], 'Tokelau': ['Tokelau Polynesian'], 'American Samoa': ['American Samoa Polynesian'], 'Papua New Guinea': ['Papua New Guinea Melanesian'], 'New Caledonia': ['New Caledonia Melanesian'], 'Niue': ['Niue Polynesian'], 'Kiribati': ['Kiribati Micronesian'], 'New Zealand': ['New Zealand Polynesian'], 'Fiji': ['Fiji Melanesian'], 'Nauru': ['Nauru Micronesian'], 'Cook Islands': ['Cook Islands Polynesian']}},
           'Ethnicity': ['Amerindian', 'Arab', 'Asian', 'Australian Aborigines', 'Austronesian', 'Berber', 'Black',
                         'Caucasoid', 'Hispanic', 'Jew', 'Kurd', 'Melanesian', 'Mestizo', 'Micronesian', 'Mixed',
                         'Mulatto', 'Oriental', 'Persian', 'Polynesian', 'Siberian']}
