import concurrent.futures
import time
import random

def load_file(fn, fn2):
    s=time.time()
    time.sleep(random.uniform(1,3))
    e=time.time()
    load_time=e-s
    print(f'{fn} loaded in {round(load_time,1)}s')
    print(fn2)
    return fn

def process_file(fn):
    s=time.time()
    time.sleep(random.uniform(1,3))
    e=time.time()
    load_time=e-s
    print(f'{fn} processed in {round(load_time,1)}s')
    return f'{fn} processed'

def main():
    # Create a ProcessPoolExecutor to run tasks in parallel with separate processes
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
        # Submit tasks and get Future objects
        download_futures = [executor.submit(load_file, i, i+10) for i in range(5)]
        
        process_futures = [executor.submit(process_file,df.result()) for df in download_futures]

        for future in concurrent.futures.as_completed(process_futures):
            print(future.result())
if __name__ == '__main__':
    main()