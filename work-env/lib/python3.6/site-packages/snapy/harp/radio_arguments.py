import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-d',
    action = 'store_true',
    default = False,
    help = 'fit differential observation'
    )
parser.add_argument('-i', '--input',
    required = True,
    help = 'template input file (required)'
    )
parser.add_argument('-o', '--output',
    default = '',
    help = 'output id'
    )
parser.add_argument('--qNH3',
    default = '2.7',
    help = 'ammonia abundance'
    )
parser.add_argument('--qH2O',
    default = '20.',
    help = 'water abundance'
    )
parser.add_argument('--sNH3',
    default = '0.5',
    help = 'ammonia standard deviation'
    )
parser.add_argument('--sT',
    default = '5.',
    help = 'temperature standard deviation'
    )
parser.add_argument('--zNH3',
    default = '20',
    help = 'ammonia correlation'
    )
parser.add_argument('--zT',
    default = '20',
    help = 'temperature correlation'
    )
parser.add_argument('--T0',
    default = '166',
    help = 'reference temperature'
    )
parser.add_argument('--grav',
    default = '10.',
    help = 'surface gravity'
    )
parser.add_argument('--Tmin',
    default = '110',
    help = 'minimum temperature'
    )
parser.add_argument('--nodes',
    default = '1',
    help = 'number of nodes to use'
    )
parser.add_argument('--plevel',
    default = '50:0.6:5',
    help = 'pressure divides'
    )
parser.add_argument('--pmin',
    default = '0.3',
    help = 'minimum pressure'
    )
parser.add_argument('--pmax',
    default = '100.',
    help = 'maximum pressure'
    )
parser.add_argument('--tem',
    default = '0',
    help = 'temperature anomaly'
    )
parser.add_argument('--nh3',
    default = '0',
    help = 'ammonia anomaly'
    )
parser.add_argument('--var',
    default = '2',
    help = 'inversion variables'
    )
parser.add_argument('--nlim',
    default = '0',
    help = 'number of MCMC steps'
    )
parser.add_argument('--nwalker',
    default = '1',
    help = 'number of mcmc walkers'
    )
parser.add_argument('--obs',
    default = 'none',
    help = 'observation file'
    )
parser.add_argument('--lat',
    default = '0.',
    help = 'latitude'
    )
parser.add_argument('--rgradt',
    default = '1.',
    help = 'radiative temperature gradient'
    )
parser.add_argument('--metallicity',
    default = '0.',
    help = 'metallicity applied to Na/KCl'
    )
parser.add_argument('--karpowicz_scale',
    default = '0.',
    help = 'scale karpowicz water opacity model'
    )
parser.add_argument('--hanley_power',
    default = '0.',
    help = 'scale hanley ammonia opacity model'
    )
args = vars(parser.parse_args())
if '.' not in args['input']:
    args['input'] += '.tmp'
