import sys
from clean_dbs import find_all_blocks, clean_block

raw_folder = sys.argv[1]
out_folder = sys.argv[2]
failed = []
for block in find_all_blocks(raw_folder):
    try: #one bad block (missing file, corrupt data) shouldnt stop the rest
        print(clean_block(block,out_folder))
    except Exception as error:
        print(f"FAILED {block}: {error}")
        failed.append(block)

print(f"done, {len(failed)} failed")
